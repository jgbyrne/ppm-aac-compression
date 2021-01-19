import sys
import ppm

def main():
    try: 
        filename = sys.argv[1]
    except IndexError:
        print("python encoder.py <filename>")
        return 1

    config = ppm.Configuration(5, 66, 32)
    frqs   = ppm.Frequencies(config)
    enc    = ppm.Encoder(config, frqs)

    # Build list of symbols to encode from raw bytes
    # 0..255 byte space is mapped into 0..63 symbol space
    # Shift symbols 64 and 65 are used to achieve this 
    symbols = []
    with open(filename, 'rb') as inf:
        while True:
            chunk = inf.read(2048)
            if not chunk: break
            for byte in chunk:
                s = int(byte)
                if s < 64:
                    symbols.append(64)
                    symbols.append(s)
                elif 128 <= s < 192:
                    symbols.append(64)
                    symbols.append(65)
                    symbols.append(s - 128)
                elif s >= 192:
                    symbols.append(65)
                    symbols.append(s - 192)
                else:
                    symbols.append(s - 64)
    l = len(symbols)

    symbols.append(config.eof_sym)
    for sym in symbols:
        enc.encode(sym)

    result = enc.conclude()
    out_l = len(result)

    print("Compressed {} bytes in {} (x{:.3f})".format(l, out_l, ((out_l) / l)))

    with open(filename.rsplit('.', 1)[0] + ".lz", "wb") as outf:
        outf.write(bytes(result))
    return 0

if __name__ == "__main__":
    sys.exit(main())

