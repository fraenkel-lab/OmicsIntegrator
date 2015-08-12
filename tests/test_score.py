# Test functions for forest.score()
# score(value, mu, musquared)

import sys
sys.path.insert(0, '../scripts')

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