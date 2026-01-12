import mpmath as mp

class DFA:
    """
    Moore machine:
      - states: finite set (use ints)
      - start: start state
      - trans: dict (state, symbol) -> next_state, where symbol in {0,1}
      - out: dict state -> output bit {0,1}
    We feed the binary digits of n (or n-1) to the automaton to get s_n.
    """
    def __init__(self, states, start, trans, out):
        self.states = states
        self.start = start
        self.trans = trans
        self.out = out

    def output_for_n(self, n: int, use_n_minus_1: bool = True) -> int:
        x = n - 1 if use_n_minus_1 else n
        # feed bits MSB -> LSB; ensure at least one bit
        bits = bin(x)[2:] or "0"
        st = self.start
        for ch in bits:
            sym = 1 if ch == "1" else 0
            st = self.trans[(st, sym)]
        return int(self.out[st])


def automaton_constant(dfa: DFA, prec_bits: int = 200, use_n_minus_1: bool = True) -> mp.mpf:
    """
    Compute C_A = sum_{n>=1} 2^{-n} log(1 + s_n/n)

    Truncation:
      each term <= 2^{-n} * log(1 + 1/n) < 2^{-n} * (1/n)
      tail < sum_{n>N} 2^{-n}/n < (1/N) * sum_{n>N} 2^{-n} = (1/N) * 2^{-N}
    Choose N so that (2^{-N}/N) < 2^{-prec_bits}  (very safe).
    """
    mp.mp.dps = int(prec_bits / 3.32192809489) + 20  # bits -> decimal digits + guard

    # choose N
    N = max(10, prec_bits)
    while mp.power(2, -N) / N > mp.power(2, -prec_bits):
        N += 10

    s = mp.mpf("0")
    for n in range(1, N + 1):
        sn = dfa.output_for_n(n, use_n_minus_1=use_n_minus_1)
        if sn:
            s += mp.power(2, -n) * mp.log1p(mp.mpf(1) / n)  # log(1+1/n)
        # else add 0
    return s, N


# ---- Example DFA: Thue–Morse (parity of popcount) ----
# Output s_n = popcount(n-1) mod 2
# Moore machine:
#   state 0 = even parity, state 1 = odd parity
#   reading a '1' toggles parity; reading '0' keeps.
states = {0, 1}
start = 0
trans = {
    (0, 0): 0, (0, 1): 1,
    (1, 0): 1, (1, 1): 0,
}
out = {0: 0, 1: 1}
tm = DFA(states, start, trans, out)

val, N_used = automaton_constant(tm, prec_bits=300, use_n_minus_1=True)
print("C_A =", val)
print("N_used =", N_used)
