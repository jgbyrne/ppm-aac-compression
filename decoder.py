import sys
import ppm

def main():
    try:
        filename = sys.argv[1]
    except IndexError:
        print("python decoder.py <filename>")
        return 1

    config = ppm.Configuration(5, 66, 32)
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

    # Decode compressed file
    # Remaps from 0..63 space back into 0..255 byte space
    dec = ppm.Decoder(config, frqs, bitstring)
    out_bytes = []
    decr = -64
    while True:
        sym = dec.decode()
        if sym == config.eof_sym:
            break
        elif sym == 64:
            decr += 64
            continue
        elif sym == 65:
            decr -= 128
            continue
        elif sym == None:
            print("Decoder Crash :-(")
            return 1

        out_bytes.append(sym - decr)
        decr = -64 

    with open(filename.split('.')[0] + "-decoded.tex", 'wb') as outf:
        outf.write(bytes(out_bytes))
    
if __name__ == "__main__":
    sys.exit(main())

