import sys
import ppm

def main():
    try:
        filename = sys.argv[1]
    except IndexError:
        print("python decoder.py <filename>")
        return 1

    config = ppm.Configuration(5, 256, 32)
    frqs   = ppm.Frequencies(config)
    enc    = ppm.Encoder(config, frqs)

    with open(filename, 'rb') as inf:
        in_bytes = []
        while True:
            chunk = inf.read(2048)
            if not chunk: break
            for byte in chunk:
                in_bytes.append(byte)
        bitstring = ppm.Bitstring.from_bytes(in_bytes)

    dec = ppm.Decoder(config, frqs, bitstring)
    out_bytes = []
    while True:
        sym = dec.decode()
        if sym == config.eof_sym:
            break
        if sym == None:
            print("Decoder Crash :-(")
            return 1
        out_bytes.append(sym)

    with open(filename.split('.')[0] + "-decoded.tex", 'wb') as outf:
        outf.write(bytes(out_bytes))
    
if __name__ == "__main__":
    sys.exit(main())

