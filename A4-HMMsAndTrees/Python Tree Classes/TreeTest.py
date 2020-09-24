"""
Python Tree Test
Stuart Bradley - 5931269
23-05-2014
"""

from Node import Node
from Tree import Tree

class TreeTest:
	def __init__(self): 
		pass

	def test(self):
		tree = Tree()
		root = Node(label="root")
		A = Node(label="a")
		B = Node(label="b")
		C = Node(label="c")
		D = Node(label="d")

		D.add_child(A)
		D.add_child(B)
		root.add_child(D)
		root.add_child(C)

		A.set_height(0.0)
		B.set_height(0.0)
		C.set_height(0.0)
		D.set_height(1.0)
		root.set_height(2.0)

		tree.set_root(root)

		if (tree.get_leaf_count() == 3):
			print("True")
		else:
			print(tree.get_leaf_count())

		if (tree.get_newick() == "((a:1.0,b:1.0)d:1.0,c:2.0)root:0.0;"):
			print("True")

s = TreeTest()
s.test()