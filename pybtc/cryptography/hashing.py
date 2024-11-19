from hashlib import sha256

def hash256(m):
    '''Performs two rounds of SHA256 on an input message'''
    return sha256(sha256(m).digest()).digest()

