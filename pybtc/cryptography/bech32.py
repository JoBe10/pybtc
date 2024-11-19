CHARSET = "qpzry9x8gf2tvdw0s3jn54khce6mua7l"

def bech32_polymod(values: list[int]) -> int:
    """
    Internal function to compute the Bech32 checksum.

    Args:
        values (list[int]): List of integer values for checksum computation.

    Returns:
        int: The computed checksum as an integer.
    """
    GEN = [0x3b6a57b2, 0x26508e6d, 0x1ea119fa, 0x3d4233dd, 0x2a1462b3]
    chk = 1
    for v in values:
        b = chk >> 25
        chk = ((chk & 0x1ffffff) << 5) ^ v
        for i in range(5):
            if (b >> i) & 1:
                chk ^= GEN[i]
    return chk

def bech32_hrp_expand(hrp: str) -> list[int]:
    """
    Expand the HRP for Bech32 checksum calculation.

    Args:
        hrp (str): The human-readable part.

    Returns:
        list[int]: Expanded HRP as a list of integers.
    """
    return [ord(x) >> 5 for x in hrp] + [0] + [ord(x) & 31 for x in hrp]

def bech32_create_checksum(hrp: str, data: list[int]) -> list[int]:
    """
    Compute the Bech32 checksum.

    Args:
        hrp (str): The human-readable part.
        data (list[int]): Data to encode.

    Returns:
        list[int]: The computed checksum as a list of integers.
    """
    values = bech32_hrp_expand(hrp) + data
    polymod = bech32_polymod(values + [0, 0, 0, 0, 0, 0]) ^ 1
    return [(polymod >> 5 * (5 - i)) & 31 for i in range(6)]

def bech32_encode(hrp: str, data: list[int]) -> str:
    """
    Encode a Bech32 string.

    Args:
        hrp (str): The human-readable part.
        data (list[int]): Data to encode.

    Returns:
        str: The encoded Bech32 string.
    """
    combined = data + bech32_create_checksum(hrp, data)
    return hrp + '1' + ''.join([CHARSET[d] for d in combined])

def convertbits(data: list[int], frombits: int, tobits: int, pad: bool = True) -> list[int]:
    """
    General power-of-2 base conversion.

    Args:
        data (list[int]): The data to convert.
        frombits (int): The number of bits of each input value.
        tobits (int): The number of bits of each output value.
        pad (bool): Whether to pad the output. Default is True.

    Returns:
        Optional[list[int]]: The converted data as a list of integers, or None if conversion fails.
    """
    acc = 0
    bits = 0
    ret = []
    maxv = (1 << tobits) - 1
    for value in data:
        if value < 0 or (value >> frombits):
            return None
        acc = (acc << frombits) | value
        bits += frombits
        while bits >= tobits:
            bits -= tobits
            ret.append((acc >> bits) & maxv)
    if pad:
        if bits:
            ret.append((acc << (tobits - bits)) & maxv)
    elif bits >= frombits or ((acc << (tobits - bits)) & maxv):
        return None
    return ret

def encode_segwit_address(hrp: str, witver: int, witprog: bytes) -> str:
    """
    Encode a SegWit address using Bech32 encoding.

    Args:
        hrp (str): The human-readable part (e.g., 'bc' for mainnet or 'tb' for testnet).
        witver (int): The witness version.
        witprog (bytes): The witness program.

    Returns:
        str: The encoded SegWit address.
    """
    data = [witver] + convertbits(list(witprog), 8, 5)
    return bech32_encode(hrp, data)
