# Test functions for forest.score()
# score(value, mu, musquared)

import os, sys

# Create the path to forest relative to the test_score.py path
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

from forest import score

class TestScore:
	
	def test_muGreaterThanZero(self):
		try:
			score(2, 0, 0)
		except ValueError:
			assert 1
		else:
			assert 0
			
		try:
			score(2, -1, 0)
		except ValueError:
			assert 1
		else:
			assert 0
		
	def test_whenValueOne(self):
		assert score(1, 1, 0) == 0
		
	def test_musquaredFalse(self):
		assert score(2, 2, 0) == -4
		
	def test_musquaredTrue(self):
		assert score(2, 2, 1) == -8
