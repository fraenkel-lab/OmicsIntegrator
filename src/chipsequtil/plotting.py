import math

from matplotlib.pyplot import hist, plot, savefig, title, show, xticks, yticks, figure, clf

from chipsequtil import get_gc_content

def plot_gc_content(sequences,bins=10,fn=None) :

    # calculate all the GC contents, sort them
    gc_contents = map(get_gc_content,sequences)
    gc_contents.sort()

    f = figure()
    points = hist(gc_contents,bins=bins)
    if fn :
        savefig(fn)
    else :
        show()
    clf()


def plot_pos_neg_peaks(pos_peaks,neg_peaks) :
    '''Plot # pos peaks/# neg peaks by p-value'''
    pass
