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
# : alphabet     - size of alphabet (excluding EOF, ESC)
# : magnitude    - number of working fixed-point places 
class Configuration:
    def __init__(self, initial_order, alphabet, magnitude):
        self.initial_order   = initial_order
        self.eof_sym         = alphabet
        self.esc_sym         = alphabet + 1

        self.normal_symbols  = alphabet + 1
        self.symbol_count    = alphabet + 2

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

    def record(self, ctx, sym):
        ofrq = self.frq[len(ctx)]

        ctx_map = ofrq.get(ctx)
        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        if sym in ctx_map: 
            ctx_map[sym] += 1
        else:
            ctx_map[sym] = 1 
            ctx_map[self.config.esc_sym] += 1 

    # Get interval for given (ctx, sym) - used for encoding
    def interval(self, order, ctx, sym, exclude):
        if order == -1:
            return (True, sym / self.config.normal_symbols, (sym + 1) / self.config.normal_symbols)

        ofrq      = self.frq[order]
        ctx_map   = ofrq.get(ctx)

        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        acc   = 0
        left  = 0
        right = 0
        for i, (s, f) in enumerate(ctx_map.items()):
            if s != self.config.esc_sym:
                if s in exclude:
                    continue
                exclude.add(s)
            if s == sym:
                left  = acc
                acc = right = acc + f
            else:
                acc += f

        if not right:
            return (False, exclude)

        return (True, left/acc, right/acc)

    # Get interval for given (ctx, pt) - used for decoding
    def query(self, order, ctx, pt, exclude):
        if order == -1:
            sym  = floor(pt * self.config.normal_symbols)
            return (sym, sym / self.config.normal_symbols, (sym + 1) / self.config.normal_symbols, exclude)

        ofrq    = self.frq[order]
        ctx_map = ofrq.get(ctx)

        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        acc   = 0
        pts = []
        for s, f in ctx_map.items():
            if s != self.config.esc_sym:
                if s in exclude:
                    continue
                exclude.add(s)
            nxt = acc + f
            pts.append((s, acc, nxt))
            acc = nxt 

        target = pt * acc
        for pt in pts:
            if target < pt[2]:
                return (pt[0], pt[1]/acc, pt[2]/acc, exclude)

def sub_ctx(ctx, order):
    if order < 1: 
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
    def _encode_symbol(self, order, ctx, sym, exclude):
        sym_int = self.frqs.interval(order, ctx, sym, exclude)
        if not sym_int[0]:
            return sym_int # new exclude

        _, low_pt, high_pt = sym_int

        span = self.high - self.low
        self.high = self.low + floor(span * high_pt)  - 1
        self.low  = self.low + floor(span * low_pt) + 1

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

        return (True,)

    # External interface to encode one symbol
    # : Automatically cascades through orders
    def encode(self, sym):
        order = top_order = len(self.ctx)

        exclude = set()
        while order >= -1:
            ctx = sub_ctx(self.ctx, order) 
            enc_res = self._encode_symbol(order, ctx, sym, exclude.copy())
            if enc_res[0]:
                break
            else:
                self._encode_symbol(order, ctx, self.config.esc_sym, exclude)
                #self.frqs.record(sub_ctx(self.ctx, order), self.config.esc_sym)

            exclude = enc_res[1]
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
    def _decode_symbol(self, order, ctx, exclude):
        span = self.high - self.low
        pt = (self.num - self.low) / span
        sym, low_pt, high_pt, exclude = self.frqs.query(order, ctx, pt, exclude) 

        self.high = self.low + floor(span * high_pt) - 1
        self.low  = self.low + floor(span * low_pt) + 1

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

        return (sym, exclude)

    # External interface to decode one symbol
    # : Automatically cascades through orders
    def decode(self):
        order = top_order = len(self.ctx)
        exclude = set()
        while order >= -1:
            ctx = sub_ctx(self.ctx, order)
            sym, exclude = self._decode_symbol(order, ctx, exclude)
            if sym != self.config.esc_sym:
                break
            order -= 1

        for o in range(top_order + 1):
            self.frqs.record(sub_ctx(self.ctx, o), sym)
   
        self.ctx.append(sym)
        if len(self.ctx) > self.config.initial_order:
            self.ctx.pop(0)

        return sym

# Just for fun :-)
if __name__ == "__main__": 
    bs=Bitstring.from_bytes([194,163,205,236,240,108,218,96,210,116,137,86,171,150,190,
                             125,9,34,224,133,1,247,255,226,95,112,122,42,135,216,224,
                             51,87,161,161,145,237,198,15,173,26,224,196,216,43,248])
    c=Configuration(5,256,32);f=Frequencies(c);d=Decoder(c,f,bs);b=[]
    while True:
        sym = d.decode()
        if sym == c.eof_sym: print("".join(bytes(b).decode("utf8"))); break
        b.append(sym)

