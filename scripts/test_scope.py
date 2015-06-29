# Test functions for forest.scope()
# scope(self, value, mu, musquared)

from forest import score

class TestScope:
 
	def test_whenValueOne(self):
		assert score(self, 1, 0, 0) == 0
		
	def test_musquaredFalse(self):
		assert score(self, 2, 2, 0) == -4
		
	def test_musquaredTrue(self):
		assert score(self, 2, 2, 1) == -8