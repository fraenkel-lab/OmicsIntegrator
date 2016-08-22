import copy, os, shlex, shutil, subprocess, sys, tempfile

# import loadGraph from scripts utilities
cur_dir = os.path.dirname(__file__)
path = os.path.abspath(os.path.join(cur_dir, '..', 'scripts'))
if not path in sys.path:
    sys.path.insert(1, path)
del path
from forest_util import loadGraph

def print_graph(g):
    '''
    Print graph g for resolving test assertion errors
    '''
    print("Printing graph:")
    for e in g.edges_iter():
      print(e)

def write_conf(fh, params):
    '''
    Write a file to use for the --conf option in forest.py
    
    INPUT:
    fh - a file-like object
    params - the values of w, b, D and optionally mu, garnetBeta, noise, r, g.
    '''
    for k,v in params.iteritems():
        if k == 'D' or k == 'processes':
          fh.write('%s = %d\n' % (k, v))
        else:
          fh.write('%s = %f\n' % (k, v))

def parse_obj(info_filename):
    '''
    Parse the objective function value from the info file.  Raises an exception
    if the objective function value is not found.
    
    INPUT:
    info_filename - the path of the info file
    
    OUTPUT:
    objective - the objective function value of the optimal Steiner forest (float)
    '''
    with open(info_filename) as info_f:
        for line in info_f:
            if line.startswith('Total objective function:'):
                parts = line.split(':')
                return float(parts[1].strip())
        
        # The objective function value was not found in the info file
        raise RuntimeError('Objective function value not found in %s' % info_filename)

def run_forest(msgsteiner, conf_params, forest_opts):
    '''
    Run Forest on the test network
    
    INPUT:
    msgsteiner - the path to msgsteiner
    conf_params - dictionary, see write_conf
    forest_opts - dictionary with long form option names for forest.py as keys

    OUTPUT:
    graph - the DiGraph object for the optimal Steiner forest
    objective - the objective function value of the optimal Steiner forest (float)
    '''
    assert msgsteiner is not None, 'Please provide path to msgsteiner using --msgpath option'
    forest_opts = copy.deepcopy(forest_opts)
    forest_opts['msgpath'] = msgsteiner

    # For reproducibility, though it is likely not needed for this
    # small test case        
    forest_opts['seed'] = 2016
    
    # Write the configuration file with the specified parameters
    conf_filename= tempfile.NamedTemporaryFile(suffix='.txt', delete=False)
    try:
        write_conf(conf_filename, conf_params)
    finally:
        conf_filename.close()
    forest_opts['conf'] = conf_filename.name
    
    default_outpath = True
    if 'outpath' in forest_opts:
      default_outpath = False
    else:
      # Create a tmp directory for output unless one is provided
      forest_opts['outpath'] = tempfile.mkdtemp()
    
    try:
        cur_dir = os.path.dirname(__file__)
        forest_opts['forest'] = os.path.join(cur_dir, '..', 'scripts', 'forest.py')
        forest_cmd = 'python {forest} --prize={prize} --edge={edge} --conf={conf} --dummyMode={dummyMode} --outpath={outpath} --msgpath={msgpath} --seed={seed}'.format(**forest_opts)
        subprocess.call(shlex.split(forest_cmd), shell=False)	

        # Test the optimal Forest to see if parameter value has the intended effect
        opt_forest = os.path.join(forest_opts['outpath'], 'result_optimalForest.sif')
        assert os.path.isfile(opt_forest), 'Forest did not generate the optimal forest file'
        graph = loadGraph(opt_forest)
        
        # Parse the objective function value of the optimal forest
        # This cannot be recovered from the forest file because it depends on w
        info_filename = os.path.join(forest_opts['outpath'], 'result_info.txt')
        assert os.path.isfile(info_filename), 'Forest did not generate the info file'
        objective = parse_obj(info_filename)
        
    finally:
        # Remove here because delete=False above
        os.remove(conf_filename.name)
        # Remove the Forest output directory and all files
        # Leave the output directory intact if it was manually set
        if default_outpath:
          shutil.rmtree(forest_opts['outpath'])
        
    return (graph, objective)
