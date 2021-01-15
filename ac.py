from math import floor

class Bitstring:
    def __init__(self):
        self.string = ""

    def push(self, val):
        if val not in (0, 1):
            raise ValueError
        self.string += str(val)
    
    def pop(self):
        if not self.string:
            return 0
        val = int(self.string[0])
        self.string = self.string[1:]
        return val

    def __repr__(self):
        return self.string

    def __str__(self):
        return self.string

    def __len__(self):
        return len(self.string)

class Frequencies:
    def __init__(self, frq):
        self.frq = frq

    def interval(self, sym):
        total = 0
        left  = 0
        right = 0
        for i, f in enumerate(self.frq):
            if i == sym:
                left = total
                right = total + f
            total += f
        return (left/total, right/total)

    def query(self, pt):
        fsum = sum(self.frq)
        point = pt * fsum
        acc = 0
        for i, f in enumerate(self.frq):
            nxt = acc + f
            if nxt > point:
                return (i, acc / fsum, nxt / fsum)
            acc = nxt

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
    config = Configuration(16)
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

def tex():
    flist = []
    with open("counts") as inf:
        for line in inf:
            flist.append(int(line.strip()))
    flist = flist[:256]
    flist.append(2)
    frqs = Frequencies(flist)
    config = Configuration(24)
    enc = Encoder(config)

    with open("in.tex") as inf:
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
    tex()

