# Test functions for forest.score()
# score(value, mu, musquared)

import sys
sys.path.insert(0, '../scripts')

from forest import score

class TestScore:
 
	def test_whenValueOne(self):
		assert score(1, 0, 0) == 0
		
	def test_musquaredFalse(self):
		assert score(2, 2, 0) == -4
		
	def test_musquaredTrue(self):
		assert score(2, 2, 1) == -8