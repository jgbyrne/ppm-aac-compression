from math import floor

# Bitstring
# | Stores encoded data
# : Backing structure is a list of Integers in range [0, 255]
class Bitstring:
    def __init__(self):
        self.bytes = []
        self.rbit = 7
        self.lbit = 7

    @classmethod
    def from_bytes(cls, blist):
        b = cls()
        b.bytes = blist
        return b

    def push(self, val):
        if val not in (0, 1):
            raise ValueError
        if not self.bytes or self.rbit == 0:
            self.bytes.append(val << 7)
            self.rbit = 7
        else:
            self.rbit -= 1
            self.bytes[-1] += (val << self.rbit)

    def pop(self):
        if not self.bytes:
            return 0
        pop_bit = int(bool(self.bytes[0] & (1 << self.lbit)))
        if self.lbit == 0:
            self.bytes.pop(0)
            self.lbit = 7
        else:
            self.lbit -= 1
        return pop_bit

    def __repr__(self):
        return "".join("{:08b}".format(b) for b in self.bytes)

    def __str__(self):
        return "".join("{:08b}".format(b) for b in self.bytes)

    def __len__(self):
        return len(self.bytes) * 8

# Configuration
# | Static paramaters for Encoding and Decoding
# : intial_order - maximum context length
# : num_symbols  - size of alphabet (excluding EOF, ESC)
# : magnitude    - number of working fixed-point places 
class Configuration:
    def __init__(self, initial_order, num_symbols, magnitude):
        self.initial_order   = initial_order
        self.num_symbols     = num_symbols
        self.eof_sym         = num_symbols
        self.esc_sym         = num_symbols + 1

        self.magnitude       = magnitude
        self.max             = 2 ** magnitude
        self.half            = self.max >> 1
        self.quarter         = self.max >> 2
        self.three_quarters  = self.half + self.quarter

        self.full_mask       = 2 ** magnitude
        self.half_mask       = 2 ** (magnitude - 1)
        self.chunk_mask      = self.full_mask - 1

# Frequencies
# | Adaptive table of (context, symbol) frequencies 
class Frequencies:
    def __init__(self, config):
        self.config = config
        self.frq    = []
        for i in range(config.initial_order + 1):
            self.frq.append({})

    def populate(self):
        for i in range(self.config.num_symbols):
            self.overwrite((), i, 1)
        self.overwrite((), self.config.eof_sym, 1)

    def overwrite(self, ctx, sym, f):
        ofrq = self.frq[len(ctx)]

        ctx_map = ofrq.get(ctx)
        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        ctx_map[sym] = f

    def record(self, ctx, sym):
        ofrq = self.frq[len(ctx)]

        ctx_map = ofrq.get(ctx)
        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        if sym in ctx_map: 
            ctx_map[sym] += 1
        else:
            ctx_map[sym] = 1

    # Get interval for given (ctx, sym) - used for encoding
    # : order is implied by context length
    def interval(self, ctx, sym):
        order     = len(ctx)
        ofrq      = self.frq[order]
        ctx_map   = ofrq.get(ctx)

        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        acc   = 0
        left  = 0
        right = 0
        for i, (s, f) in enumerate(ctx_map.items()):
            if s == sym:
                left  = acc
                acc = right = acc + f
            else:
                acc += f

        if not right:
            return None

        return (left/acc, right/acc)

    # Get interval for given (ctx, pt) - used for decoding
    # : order is implied by context length
    def query(self, ctx, pt):
        order   = len(ctx)
        ofrq    = self.frq[order]
        ctx_map = ofrq.get(ctx)

        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        total  = sum(ctx_map.values())
        target = pt * total

        acc   = 0
        left  = 0
        right = 0
        sym   = None
        for i, (sym, f) in enumerate(ctx_map.items()):
            nxt = acc + f
            if nxt > target:
                left  = acc
                right = nxt
                break
            acc = nxt

        return (sym, left/total, right/total)

def sub_ctx(ctx, order):
    if not order: 
        return ()
    else:
        return tuple(ctx[-order:])

# Encoder
# | Produces a PPM encoding from a sequence of symbols
class Encoder:
    def __init__(self, config, frqs):
        self.config   = config
        self.frqs     = frqs

        self.low      = 0
        self.high     = config.max - 1

        self.ctx      = []

        self.straddle = 0
        self.bstr     = Bitstring()

    # Internal symbol encoding routine
    # : If this fails, sym_esc is encoded and the order is dropped
    def _encode_symbol(self, ctx, sym):
        sym_int = self.frqs.interval(ctx, sym)
        if sym_int is None:
            return False

        low_pt, high_pt = sym_int

        span = self.high - self.low
        self.high = self.low + floor(span * high_pt)  - 1
        self.low  = self.low + floor(span * low_pt) + 1
        
        assert (self.high - self.low) > 0

        while True:
            if not self.high & self.config.half_mask:
                self.bstr.push(0)
                for _ in range(self.straddle):
                    self.bstr.push(1)
                self.straddle = 0
                self.low  = self.low << 1
                self.high = (self.high << 1) + 1

            elif self.low & self.config.half:
                self.bstr.push(1)
                for _ in range(self.straddle):
                    self.bstr.push(0)
                self.straddle = 0
                self.low  = (self.low << 1) ^ self.config.full_mask
                self.high = ((self.high << 1) ^ self.config.full_mask) + 1

            elif self.config.quarter <= self.low and self.high < self.config.three_quarters:
                self.straddle += 1
                self.high = ((self.high << 1) ^ (self.config.full_mask | self.config.half_mask)) + 1
                self.low  = (self.low << 1) ^ self.config.half_mask

            else:
                break

        return True

    # External interface to encode one symbol
    # : Automatically cascades through orders
    # : N.B. Order 0 encoding is assumed to never fail
    def encode(self, sym):
        order = top_order = len(self.ctx)

        while order >= 0:
            ctx = () if not order else tuple(self.ctx[-order:])
            if self._encode_symbol(ctx, sym):
                break
            else:
                self._encode_symbol(ctx, self.config.esc_sym)
                self.frqs.record(sub_ctx(self.ctx, order), self.config.esc_sym)
            order -= 1

        for o in range(top_order + 1):
            self.frqs.record(sub_ctx(self.ctx, o), sym)

        self.ctx.append(sym)
        if len(self.ctx) > self.config.initial_order:
            self.ctx.pop(0)

    # External interface to finalise encoding and produce byte list
    # : Final encoding is midpoint between high and low
    def conclude(self):
        mid = self.low + ((self.high - self.low) // 2)
        if self.straddle:
            if mid & self.config.half_mask:
                self.bstr.push(1)
                for _ in range(self.straddle):
                    self.bstr.push(0)
                mid = (mid << 1) ^ self.config.full_mask
            else:
                self.bstr.push(0)
                for _ in range(self.straddle):
                    self.bstr.push(1)
                mid = mid << 1

        while mid:
            if mid & self.config.half_mask:
                self.bstr.push(1)
                mid = (mid << 1) ^ self.config.full_mask
            else:
                self.bstr.push(0)
                mid = (mid << 1)
        return self.bstr.bytes

# Encoder
# | Produces a sequence of symbols from a PPM encoding
class Decoder:
    def __init__(self, config, frqs, bstr):
        self.config   = config
        self.frqs     = frqs
        
        self.low      = 0
        self.high     = config.max - 1

        self.ctx      = []

        self.straddle = 0
        self.bstr     = bstr
        self.num      = 0
        
        for _ in range(self.config.magnitude):
            self.shift()

    # Reproduce a Lower or Upper normalisation
    def shift(self):
        nxt_bit = self.bstr.pop()
        self.num <<= 1
        self.num &= self.config.chunk_mask
        self.num += nxt_bit

    # Reproduce a Middle normalisation
    def inflate(self):
        nxt_bit = self.bstr.pop()
        self.num <<= 1
        side = self.num & self.config.full_mask
        if side:
            self.num = (self.num - self.config.max) + self.config.half + nxt_bit
        else:
            self.num = (self.num - self.config.half) + nxt_bit

    # Internal symbol decoding routine
    def _decode_symbol(self, ctx):
        span = self.high - self.low
        pt = (self.num - self.low) / span
        sym, low_pt, high_pt = self.frqs.query(ctx, pt) 

        self.high = self.low + floor(span * high_pt) - 1
        self.low  = self.low + floor(span * low_pt) + 1
        
        assert (self.high - self.low) > 0

        while True:
            if not self.high & self.config.half_mask:
                self.shift()
                self.low  = self.low << 1
                self.high = (self.high << 1) + 1

            elif self.low & self.config.half:
                self.shift()
                self.low  = (self.low << 1) ^ self.config.full_mask
                self.high = ((self.high << 1) ^ self.config.full_mask) + 1

            elif self.config.quarter <= self.low and self.high < self.config.three_quarters:
                self.high = ((self.high << 1) ^ (self.config.full_mask | self.config.half_mask)) + 1
                self.low  = (self.low << 1) ^ self.config.half_mask
                self.inflate()

            else:
                break

        return sym

    # External interface to decode one symbol
    # : Automatically cascades through orders
    # : N.B. Order 0 decoding is assumed to never fail
    def decode(self):
        order = top_order = len(self.ctx)

        while order >= 0:
            ctx = sub_ctx(self.ctx, order)
            sym = self._decode_symbol(ctx)
            if sym != self.config.esc_sym:
                break
            else:
                self.frqs.record(sub_ctx(self.ctx, order), self.config.esc_sym)
            order -= 1

        for o in range(top_order + 1):
            self.frqs.record(sub_ctx(self.ctx, o), sym)
   
        self.ctx.append(sym)
        for i in range(len(self.ctx) - self.config.initial_order):
            self.ctx.pop(0)

        return sym

# Just for fun :-)
if __name__ == "__main__": 
    bs=Bitstring.from_bytes([195,19,125,179,112,236,70,85,217,133,85,36,228,173,
                             93,95,249,219,92,237,37,126,215,42,228,29,46,31,188,
                             176,178,255,143,46,104,49,99,63,122,222,187,105,110,
                             27,209,109,104,120,181,252])
    c=Configuration(4,256,32);f=Frequencies(c);d=Decoder(c,f,bs);f.populate();b=[]
    while True:
        sym = d.decode()
        if sym == c.eof_sym: print("".join(bytes(b).decode("utf8"))); break
        b.append(sym)

