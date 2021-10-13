
def float_to_ieee(num):
    import struct
    bits, = struct.unpack('!I', struct.pack('!f', num))
    return "{:032b}".format(bits)

def ieee_to_float(N): # ieee-745 bits (max 32 bit)
    a = int(N[0])        # sign,     1 bit
    b = int(N[1:9],2)    # exponent, 8 bits
    c = int("1"+N[9:], 2)# fraction, len(N)-9 bits

    return (-1)**a * c /( 1<<( len(N)-9 - (b-127) ))
