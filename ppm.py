from math import floor
import pprofile
import copy
import pprint

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

class Frequencies:
    def __init__(self, config):
        self.config = config
        self.frq    = []
        for i in range(config.initial_order + 1):
            self.frq.append({})

    def overwrite(self, ctx, sym, f):
        ofrq = self.frq[len(ctx)]

        if (ctx_map := ofrq.get(ctx)) is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        ctx_map[sym] = f

    def record(self, ctx, sym, debug=False):
        if debug:
            print("Recording {} with context {}".format(sym, ctx))
        ofrq = self.frq[len(ctx)]

        if (ctx_map := ofrq.get(ctx)) is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        if sym in ctx_map:
            ctx_map[sym] += 1
        else:
            ctx_map[sym] = 1

    def interval(self, ctx, sym):
        order     = len(ctx)
        ofrq      = self.frq[order]
        ctx_map   = ofrq.get(ctx)

        if ctx_map is None:
            ctx_map = ofrq[ctx] = {self.config.esc_sym: 1}

        acc   = 0
        left  = 0
        right = 0
        for i, sym_frq in enumerate(ctx_map.items()):
            if sym_frq[0] == sym:
                left  = acc
                acc = right = acc + sym_frq[1]
            else:
                acc += sym_frq[1]

        if not right:
            return None

        return (left/acc, right/acc)

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
        for i, sym_frq in enumerate(ctx_map.items()):
            nxt = acc + sym_frq[1]
            if nxt > target:
                sym   = sym_frq[0]
                left  = acc
                right = nxt
                break
            acc = nxt

        return (sym, left/total, right/total)

class Encoder:
    def __init__(self, config, frqs):
        self.config   = config
        self.frqs     = frqs

        self.low      = 0
        self.high     = config.max - 1

        self.ctx      = []

        self.straddle = 0
        self.bstr     = Bitstring()

    def encode_symbol(self, ctx, sym, debug=False):
        #print("Encode({})".format(len(ctx)), sym, )
        sym_int = self.frqs.interval(ctx, sym)
        #print(sym_int)
        if sym_int is None:
            return False

        low_pt, high_pt = sym_int
        span = self.high - self.low
        self.high = self.low + floor(span * high_pt)  - 1
        self.low  = self.low + floor(span * low_pt)

        if debug:
            print("Encoder({})\t\t{}\t\t{}\t{}".format(len(ctx), sym, self.low, self.high))

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

    def encode(self, sym, debug=False):
        #print(self.ctx)
        order = len(self.ctx)
        if debug: print("Encode", sym, order)
        done = False
        while order >= 0:
            ctx = () if not order else tuple(self.ctx[-order:])
            if not done:
                success = self.encode_symbol(ctx, sym, debug=debug)
                if success:
                    done = True
                else:
                    self.encode_symbol(ctx, self.config.esc_sym)
                    self.frqs.record(ctx, self.config.esc_sym)
              
            self.frqs.record(ctx, sym, debug=debug)
            order -= 1

        #print("Encoded: ", self.frqs.frq[1:])
        #print(self.frqs.frq[0])

        self.ctx.append(sym)
        for i in range(len(self.ctx) - self.config.initial_order):
            self.ctx.pop(0)

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
        return self.bstr

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

    def shift(self):
        nxt_bit = self.bstr.pop()
        self.num <<= 1
        self.num &= self.config.chunk_mask
        self.num += nxt_bit

    def inflate(self):
        nxt_bit = self.bstr.pop()
        self.num <<= 1
        side = self.num & self.config.full_mask
        if side:
            self.num = (self.num - self.config.max) + self.config.half + nxt_bit
        else:
            self.num = (self.num - self.config.half) + nxt_bit

    def decode_symbol(self, ctx, debug=False):
        #print("Decode({})".format(len(ctx)), self.num, self.frqs.frq[1:])
        span = self.high - self.low
        pt = (self.num - self.low) / span
        sym, low_pt, high_pt = self.frqs.query(ctx, pt) 
        #print(low_pt, high_pt)

        self.high = self.low + floor(span * high_pt) - 1
        self.low  = self.low + floor(span * low_pt)

        if debug:
            print("Decoder({})\t{}\t\t{}\t{}".format(len(ctx), self.num, self.low, self.high))

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

    def sub_ctx(self, o):
        if o == 0:
            return ()
        elif 0 < o <= len(self.ctx):
            return tuple(self.ctx[-o:])
        else:
            print("???", o)
            raise ValueError

    def decode(self, debug=False):
        if debug: print("Decode", self.num)
        order = top_order = len(self.ctx)

        while order >= 0:
            ctx = self.sub_ctx(order)
            sym = self.decode_symbol(ctx, debug=debug)
            if sym != self.config.esc_sym:
                for o in range(top_order, -1, -1):
                    self.frqs.record(self.sub_ctx(o), sym, debug=debug)
                break
            else:
                self.frqs.record(ctx, self.config.esc_sym)
            order -= 1
        
        self.ctx.append(sym)
        for i in range(len(self.ctx) - self.config.initial_order):
            self.ctx.pop(0)

        return sym

def tex(filename):
    flist = []
    with open("counts") as inf:
        for line in inf:
            flist.append(int(line.strip()))
    flist = flist[:256]
    flist.append(2)

    config = Configuration(4, 256, 64)
    frqs = Frequencies(config)

    for i, f in enumerate(flist):
        frqs.overwrite((), i, f)

    check = 42 

    enc = Encoder(config, frqs)
    frc = None
    in_length = 0
    with open(filename, "rb") as inf:
        while (chunk := inf.read(2048)):
            for byte in chunk:
                enc.encode(int(byte), debug = False)
                in_length += 1
        enc.encode(config.eof_sym)

    #for byte in bytes("Hello! How are you today?", "utf8"):
    #    enc.encode(int(byte))
    #enc.encode(config.eof_sym)

    result = enc.conclude()

    out_length = len(result)
    print("Encoded {} bytes in {} bits ({:.3f})".format(in_length,
                                                        out_length,
                                                        ((out_length / 8) / in_length)))

    with open("encoded.lz", "wb") as outf:
        outf.write(bytes(result.bytes))

    del result
    del frqs

    frqs = Frequencies(config)

    for i, f in enumerate(flist):
        frqs.overwrite((), i, f)

    with open("encoded.lz", "rb") as inf:
        blist = []
        while (chunk := inf.read(2048)):
            for byte in chunk:
                blist.append(byte)
        
        bstr = Bitstring.from_bytes(blist)

    dec = Decoder(config, frqs, bstr)
    symbols = []
    i = 1
    while True:
        if i % 131072 == 0:
            print(i)
        sym = dec.decode()
        if sym == config.eof_sym or sym is None:
            break
        symbols.append(sym)
        i += 1

    with open("out.tex", "wb") as outf:
        outf.write(bytes(symbols))

if __name__ == "__main__":
    tex("test.tex")

