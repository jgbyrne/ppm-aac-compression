from math import floor
import pprofile

class Bitstring:
    def __init__(self):
        self.bytes = []
        self.rbit = 7
        self.lbit = 7

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


class Frequencies:
    def __init__(self, frq):
        self.frq = frq
        self.num = len(frq)
        self.intervals = {}
        self.lefts = []
        self.recalc()

    def recalc(self):
        self.total = sum(self.frq)
        acc = 0
        nxt = 0
        for i, f in enumerate(self.frq):
            nxt += f
            self.intervals[i] = (acc/self.total, nxt/self.total)
            self.lefts.append(acc)
            acc = nxt


    def interval(self, sym):
        return self.intervals[sym]

    def query(self, pt):
        fsum = sum(self.frq)
        point = pt * fsum
        acc = 0

        ptr = self.num >> 1
        jmp = ptr >> 1
        while jmp > 4:
            if self.lefts[ptr] > point:
                ptr -= jmp
            else:
                ptr += jmp
            jmp >>= 1

        if self.lefts[ptr] > point:
            while self.lefts[ptr] > point:
                ptr -= 1
            iv = self.intervals[ptr]
            return (ptr, iv[0], iv[1])
        else:
            while self.lefts[ptr] <= point:
                last = ptr
                ptr += 1
                if ptr == self.num:
                    break
            iv = self.intervals[last]
            return (last, iv[0], iv[1])

class Configuration:
    def __init__(self, order):
        self.order = order
        self.max = 2 ** order

class Decoder:
    def __init__(self, config, bstr):
        self.config   = config
        
        self.low      = 0
        self.high     = config.max - 1

        self.half     = config.max >> 1
        self.quarter  = config.max >> 2
        self.three_quarters = self.half + self.quarter 

        self.full_mask  = 2 ** self.config.order
        self.half_mask  = 2 ** (self.config.order - 1)
        self.chunk_mask = self.full_mask - 1

        self.straddle = 0
        self.bstr     = bstr
        self.num      = 0

        for _ in range(self.config.order):
            self.shift()

    def shift(self):
        nxt_bit = self.bstr.pop()
        self.num <<= 1
        self.num &= self.chunk_mask
        self.num += nxt_bit

    def inflate(self):
        nxt_bit = self.bstr.pop()
        self.num <<= 1
        side = self.num & self.full_mask
        if side:
            self.num = (self.num - self.config.max) + self.half + nxt_bit
        else:
            self.num = (self.num - self.half) + nxt_bit

    def decode(self, frqs):
        span = self.high - self.low
        pt = (self.num - self.low) / span
        #print(self.low, self.num, self.high, pt)
        sym, low_pt, high_pt = frqs.query(pt)

        self.high = self.low + floor(span * high_pt) - 1
        self.low  = self.low + floor(span * low_pt)

        #print("{:0{w}b} => {}\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.num, sym, self.low, self.high, self.straddle, self.bstr, w=self.config.order))
        while True:
            if not self.high & self.half_mask:
                self.shift()

                self.low  = self.low << 1
                self.high = (self.high << 1) + 1

                #print("Zoom Lower\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.low, self.high, self.straddle, self.bstr, w=self.config.order))

            elif self.low & self.half:
                self.shift()

                self.low  = (self.low << 1) ^ self.full_mask
                self.high = ((self.high << 1) ^ self.full_mask) + 1

                #print("Zoom Upper\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.low, self.high, self.straddle, self.bstr, w=self.config.order))

            elif self.quarter <= self.low and self.high < self.three_quarters:
                self.high = ((self.high << 1) ^ (self.full_mask | self.half_mask)) + 1
                self.low  = (self.low << 1) ^ self.half_mask

                self.inflate()

                #print("Zoom Mid\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.low, self.high, self.straddle, self.bstr, w=self.config.order))

            else:
                break
        return sym

class Encoder:
    def __init__(self, config):
        self.config   = config

        self.low      = 0
        self.high     = config.max - 1

        self.half     = config.max >> 1
        self.quarter  = config.max >> 2
        self.three_quarters = self.half + self.quarter

        self.full_mask = 2 ** self.config.order
        self.half_mask  = 2 ** (self.config.order - 1)

        self.straddle = 0
        self.bstr     = Bitstring()

    def encode(self, sym, frqs):
        low_pt, high_pt = frqs.interval(sym)
        span = self.high - self.low
        self.high = self.low + floor(span * high_pt)  - 1
        self.low  = self.low + floor(span * low_pt)

        #print("Encode {}\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(sym, self.low, self.high, self.straddle, self.bstr, w=self.config.order))

        while True:
            if not self.high & self.half_mask:
                self.bstr.push(0)
                for _ in range(self.straddle):
                    self.bstr.push(1)
                self.straddle = 0

                self.low  = self.low << 1
                self.high = (self.high << 1) + 1

                #print("Zoom Lower\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.low, self.high, self.straddle, self.bstr, w=self.config.order))

            elif self.low & self.half:
                self.bstr.push(1)
                for _ in range(self.straddle):
                    self.bstr.push(0)
                self.straddle = 0

                self.low  = (self.low << 1) ^ self.full_mask
                self.high = ((self.high << 1) ^ self.full_mask) + 1

                #print("Zoom Upper\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.low, self.high, self.straddle, self.bstr, w=self.config.order))

            elif self.quarter <= self.low and self.high < self.three_quarters:
                self.straddle += 1
                self.high = ((self.high << 1) ^ (self.full_mask | self.half_mask)) + 1
                self.low  = (self.low << 1) ^ self.half_mask

                #print("Zoom Mid\t{:0{w}b}\t{:0{w}b}\t{}\t{}".format(self.low, self.high, self.straddle, self.bstr, w=self.config.order))

            else:
                break

    def conclude(self):
        mid = self.low + ((self.high - self.low) // 2)
        if self.straddle:
            if mid & self.half_mask:
                self.bstr.push(1)
                for _ in range(self.straddle):
                    self.bstr.push(0)
                mid = (mid << 1) ^ self.full_mask
            else:
                self.bstr.push(0)
                for _ in range(self.straddle):
                    self.bstr.push(1)
                mid = mid << 1

        while mid:
            if mid & self.half_mask:
                self.bstr.push(1)
                mid = (mid << 1) ^ self.full_mask
            else:
                self.bstr.push(0)
                mid = (mid << 1)
        return self.bstr

def test():
    frqs = Frequencies([82, 15, 28, 43, 130, 22, 2, 61, 7, 1, 1, 4, 24, 67, 75, 19, 1, 6, 63, 91, 28, 1, 24, 1, 2, 1, 1])
    config = Configuration(24)
    enc = Encoder(config)

    intxt = "lookuponmyworksyemightyanddespair"
    print(intxt)

    inseq = [ord(c) - 97 for c in intxt] + [26]

    print("In: ", inseq)
    for sym in inseq:
        enc.encode(sym, frqs)

    result = enc.conclude()
    print(result)

    dec = Decoder(config, result)
    symbols = []
    while True:
        symbols.append(dec.decode(frqs))
        if symbols[-1] == 26:
            break

    print("Out:", symbols)
    print("".join([chr(s + 97) for s in symbols[:-1]]))

def tex(filename):
    flist = []
    with open("counts") as inf:
        for line in inf:
            flist.append(int(line.strip()))
    flist = flist[:256]
    flist.append(2)
    frqs = Frequencies(flist)
    config = Configuration(24)
    enc = Encoder(config)

    with open(filename) as inf:
        intxt = inf.read()

    inseq = [ord(c) for c in intxt] + [256]
    for i, sym in enumerate(inseq):
        enc.encode(sym, frqs)
        if i % 10000 == 0:
            print(i / len(inseq))

    result = enc.conclude()
    with open("out", "w") as outf:
        print(result, file=outf)

    in_length = len(intxt)
    out_length = len(result)
    print("Encoded {} bytes in {} bits ({:.3f})".format(in_length, out_length, ((out_length / 8) / in_length)))

    dec = Decoder(config, result)
    symbols = []
    while True:
        sym = dec.decode(frqs)
        if sym == 256:
            break
        symbols.append(sym)

    print("".join([chr(s) for s in symbols]) == intxt)

if __name__ == "__main__":
    tex("in.tex")
