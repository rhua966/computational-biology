import random
import numpy as np

# random.seed(0)


class HMM:

    STATE = ["H", "S", "T"]
    SYMBOL = ["B", "I", "N"]

    def __init__(self):

        # 3 by 3 initial np.array
        self.transition_matrix = np.reshape([15 / 16, 3 / 160, 7 / 160,
                                            1 / 20, 7 / 8, 3 / 40,
                                            1 / 12, 1 / 23, 5 / 6], (3, 3))

        self.emission_matrix = np.reshape([0.3, 0.6, 0.1,
                                           0.55, 0.15, 0.3,
                                           0.1, 0.2, 0.7], (3, 3))

        self.initial_state = int(random.uniform(0, 1) / (1/3))  # Initial state chosen with 1/3 possibility each
        self.state_sequence = ""
        self.symbol_sequence = ""

    def simulate(self, length):
        """Simulate state and symbol sequence of arbitrary length from HMM."""

        count = 0  # Length count

        # Initial condition
        state = self.initial_state
        symbol = None

        # Final sequence
        state_sequence = ""
        symbol_sequence = ""

        while count < length:

            # Generate random number to simulate the probability
            emission_probability = random.uniform(0, 1)
            transition_probability = random.uniform(0, 1)

            # Generate symbols
            if 0 < emission_probability < self.emission_matrix[state][0]:
                symbol = 0
            elif sum(self.emission_matrix[state][:1]) < emission_probability < sum(self.emission_matrix[state][:2]):
                symbol = 1
            elif sum(self.emission_matrix[state][:2]) < emission_probability < 1:
                symbol = 2
            symbol_sequence += HMM.SYMBOL[symbol]

            # Generate states
            if 0 < transition_probability < self.transition_matrix[state][0]:
                state = 0
            elif sum(self.transition_matrix[state][:1]) < transition_probability < sum(self.transition_matrix[state][:2]):
                state = 1
            elif sum(self.transition_matrix[state][:2]) < transition_probability < 1:
                state = 2
            state_sequence += HMM.STATE[state]

            count += 1

        # Set the sequence and also return
        self.state_sequence = state_sequence
        self.symbol_sequence = symbol_sequence
        return state_sequence, symbol_sequence

    def joint_probability(self, x=None, pi=None):

        if x is None: x = self.symbol_sequence
        if pi is None: pi = self.state_sequence

        state_seq = pi
        symbol_seq = x
        if len(state_seq) != len(symbol_seq): print("ERROR: Dimension Doesn't Match")

        # Initial probability
        probability = np.log((1/3) * (self.emission_matrix[HMM.STATE.index(state_seq[0])][HMM.SYMBOL.index(symbol_seq[0])]))

        for i in range(1, len(state_seq)):

            curr_state = HMM.STATE.index(state_seq[i])
            prev_state = HMM.STATE.index(state_seq[i-1])
            curr_symbol = HMM.SYMBOL.index(symbol_seq[i])

            # Emission and transmission probability
            e = self.emission_matrix[curr_state][curr_symbol]
            t = self.transition_matrix[prev_state][curr_state]
            probability += np.log(e * t)
            # print(f"i:{i}, transition from {prev_state} to {curr_state}: {t}, emission from {curr_state} to {curr_symbol}: {e}, joint probability: {probability}")

        return probability

    def forward(self, x=None):
        """Take the symbol sequence x as parameter and calculate probability using forward algorithm"""

        def log_sum(array):
            """Find the log of a sum of multiple logged number"""
            z = 0
            x0 = array[0]
            for num in array:
                z += np.exp(num - x0)
            return x0 + np.log(z)

        if x is None: x = self.symbol_sequence

        symbols = x

        # Shape = (n_states, n_symbols)
        n_rows = len(HMM.STATE)
        n_cols = len(symbols) + 1

        # Initial F matrix
        F = np.zeros((n_rows, n_cols), dtype=np.float64)
        F[:, 0] = [np.NINF, np.NINF, np.NINF]

        for i in range(1, n_cols):

            curr_symbol = HMM.SYMBOL.index(symbols[i-1])

            for j in range(n_rows):

                curr_state = j

                log_emission = np.log(self.emission_matrix[curr_state][curr_symbol])

                # First column
                if i == 1:
                    F[j, i] = log_emission + np.log(1/3)
                    continue

                log_transmission = np.log(self.transition_matrix[:, curr_state])
                logsum = log_sum(F[:, i-1] + log_transmission)

                F[j, i] = log_emission + logsum
        # print(F)
        return log_sum(F[:, -1])  # Return the logsum of the last column

    def set_state_sequence(self, states):
        self.state_sequence = states

    def set_symbol_sequence(self, symbols):
        self.symbol_sequence = symbols


HMM = HMM()
HMM.simulate(100)

# Calculate the joint probability using the 150 length sequence from Q1b
print(HMM.joint_probability())
# Calculate the joint probability using the given sequence
st = ''.join("S,S,H,H,H,T,T,S,S,S,H,H,H,H,H,H,S,S,S,S,S,S".split(','))
sy = ''.join("B,I,N,B,N,I,N,B,N,I,N,B,I,N,B,I,I,N,B,B,N,B".split(','))
print(HMM.joint_probability(sy, st))
joint_probability = HMM.joint_probability()
P_X = HMM.forward()
print(joint_probability, P_X)
