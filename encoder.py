import sys
import ppm

def main():
    try: 
        filename = sys.argv[1]
    except IndexError:
        print("python encoder.py <filename>")
        return 1

    config = ppm.Configuration(4, 256, 32)
    frqs   = ppm.Frequencies(config)
    enc    = ppm.Encoder(config, frqs)
    frqs.populate()
    
    l = 0
    with open(filename, 'rb') as inf:
        while True:
            chunk = inf.read(2048)
            if not chunk: break

            for byte in chunk:
                enc.encode(int(byte)); l += 1

        enc.encode(config.eof_sym)

    result = enc.conclude()
    out_l = len(result)

    print("Compressed {} bytes (x{:.3f})".format(l, ((out_l) / l)))

    with open(filename.split('.')[0] + ".lz", "wb") as outf:
        outf.write(bytes(result))
    return 0

if __name__ == "__main__":
    sys.exit(main())

