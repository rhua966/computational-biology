"""
Python Tree Test
Stuart Bradley - 5931269
23-05-2014
"""
from mutate import mutate
from tree69 import *
import random
import numpy as np


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
		print(tree.get_newick())

		if tree.get_leaf_count() == 3:
			print("True")
		else:
			print(tree.get_leaf_count())

		if tree.get_newick() == "((a:1.0,b:1.0)d:1.0,c:2.0)root:0.0;":
			print("True")


# s = TreeTest()
# s.test()


class Tree369:

	DNA = ["A", "C", "G", "T"]

	@staticmethod
	def simulate(n, lambd):
		"""
		A method to simulate trees according to the Yule model
		:argument
			n: number of leaves
			lambd: branching parameter lambda
		:return
			yule_tree: tree generated according to Yule Model
		"""
		k = n
		t = 0
		nodes = [Node(label=f"{x+1}") for x in range(n)]  # n initial nodes
		for node in nodes:
			node.set_height(0)
		available_nodes_index = [_ for _ in range(n)]  # Storing the index of available nodes
		# print(nodes, available_nodes_index)

		while k > 1:

			# Choose a exponential distributed number tk with rate k*lambd
			tk = Tree369.rand_exp(k*lambd)
			t += tk

			# Generate new node m at height t with child i and j
			m = Node(label=f"{2*n-k+1}")
			m.set_height(t)

			i = random.randint(0, len(available_nodes_index)-1)
			m.add_child(nodes[available_nodes_index[i]])
			del available_nodes_index[i]

			j = random.randint(0, len(available_nodes_index)-1)
			m.add_child(nodes[available_nodes_index[j]])
			del available_nodes_index[j]

			# Add node m to the available node list
			nodes.append(m)
			available_nodes_index.append(len(nodes) - 1)

			k -= 1

		# Set the root to be the remaining one node
		return Tree(nodes[available_nodes_index[0]])

	@staticmethod
	def jukes_cantor_model(t, L, mu):
		"""
		:argument
			t: tree with n leaves nodes
			L: sequence length L
			mu: mutation rate
		:return
			t: tree generated according to JC model
			seq_list: a dictionary with key being the label and value being the sequence
		"""
		seq_list = []

		def traverse(node):
			"""DFS visit nodes from tree"""

			if node.is_root():
				sequence = Tree369.rand_seq(L)
				node.set_sequence(sequence)
			else:
				branch_length = node.get_parent().get_height() - node.get_height()
				parent_sequence = node.get_parent().get_sequence()
				mutated_sequence = mutate(parent_sequence, branch_length, mu)
				node.set_sequence(mutated_sequence)

			seq_list.append((int(node.get_label()), node.get_sequence()))
			if node.is_leaf():
				return

			traverse(node.get_children()[0])
			traverse(node.get_children()[1])

		traverse(t.get_root())

		return seq_list

	@staticmethod
	def calculate_JC_distance_matrix(seq_list):
		""" Calculate the Jukes Cantor distance matrix from seq_list
		:argument
			seq_list:  dictionary of sequences to calculate distance
		:return:
			d: distance matrix with i, j entry being the distance between nodes with label (i + 1) and (j + 1)
		"""
		n_seqs = len(seq_list)  # Number of sequences
		L = len(seq_list[0][1])  # Length of sequence
		d = np.zeros((n_seqs, n_seqs), dtype=np.float64)

		for i in range(n_seqs):
			for j in range(n_seqs):
				# if i == j: continue

				# Number of different sites between two sequence
				n_diff = sum(i_char != j_char for i_char, j_char in zip(seq_list[i][1], seq_list[j][1]))
				f_ij = min((n_diff/L), 0.75-(1/L))  # Fraction of differing sites

				x = seq_list[i][0] - 1
				y = seq_list[j][0] - 1
				d[x, y] = (-3/4) * np.log(1 - (4 * f_ij / 3))

		return d

	@staticmethod
	def rand_exp(rate):
		"""
		Takes a rate parameter rate and produce an exponentially
		distributed random variable with parameter rate.
		"""
		# random.seed(0)

		# Draw a uniform distributed rv between 0 and 1
		u = random.uniform(0, 1)

		# Inverse of CDF of exponential distribution
		F_1 = lambda x: (- np.log(1 - x)) / rate

		# Transform the uniform number into exponential distribution and return
		return F_1(u)

	@staticmethod
	def rand_seq(L=10):
		"""Generate random DNA sequence of length L.
		:argument
			L: Length of sequence. Default to 10
		:return:
			sequence: np.array with entries in {0, 1, 2, 3}
		"""

		sequence = []  # init sequence

		while len(sequence) < L:
			# Generate random num in the range of (0, len(alphabet) - 1) as index,
			# append the corresponding character to the seq.
			sequence.append(random.randint(0, len(Tree369.DNA) - 1))

		return np.array(sequence)

	@staticmethod
	def convert_to_string(seq):
		"""Convert a numpy array representation of sequence to string"""
		s = ""
		for i in seq:
			s += Tree369.DNA[i]

		return s


t = Tree369.simulate(5, 0.5)

seq = Tree369.jukes_cantor_model(t, 200, 0.5)
plot_tree(t)
for label, sequence in seq:
	print(label, Tree369.convert_to_string(sequence))
distance_matrix = Tree369.calculate_JC_distance_matrix(seq)
print(distance_matrix[:5, :5])
t0 = compute_upgma_tree(distance_matrix[:5, :5])
plot_tree(t0)