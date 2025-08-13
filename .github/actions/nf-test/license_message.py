#!/usr/bin/env python3

#########################################
# Author: [DonFreed](https://github.com/DonFreed)
# File: license_message.py
# Source: https://github.com/DonFreed/docker-actions-test/blob/main/.github/scripts/license_message.py
# Source+commit: https://github.com/DonFreed/docker-actions-test/blob/aa1051a9f53b3a1e801953748d062cad74dca9a9/.github/scripts/license_message.py
# Download Date: 2023-07-04, commit: aa1051a
# This source code is licensed under the BSD 2-Clause license
#########################################

"""
Functions for generating and sending license messages
"""

# Modified from - https://stackoverflow.com/a/59835994

import argparse
import base64
import calendar
import re
import secrets
import sys

from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from datetime import datetime as dt

MESSAGE_TIMEOUT = 60 * 60 * 24  # Messages are valid for 1 day
NONCE_BYTES = 12


class DecryptionTimeout(Exception):
    # Decrypting a message that is too old
    pass


def generate_key():
    key = secrets.token_bytes(32)
    return key


def handle_generate_key(args):
    key = generate_key()
    key_b64 = base64.b64encode(key)
    print(key_b64.decode("utf-8"), file=args.outfile)


def encrypt_message(key, message):
    nonce = secrets.token_bytes(NONCE_BYTES)
    timestamp = calendar.timegm(dt.now().utctimetuple())
    data = timestamp.to_bytes(10, byteorder="big") + b"__" + message
    ciphertext = nonce + AESGCM(key).encrypt(nonce, data, b"")
    return ciphertext


def handle_encrypt_message(args):
    key = base64.b64decode(args.key.encode("utf-8"))
    message = args.message.encode("utf-8")
    ciphertext = encrypt_message(key, message)
    ciphertext_b64 = base64.b64encode(ciphertext)
    print(ciphertext_b64.decode("utf-8"), file=args.outfile)


def decrypt_message(key, ciphertext, timeout=MESSAGE_TIMEOUT):
    nonce, ciphertext = ciphertext[:NONCE_BYTES], ciphertext[NONCE_BYTES:]
    message = AESGCM(key).decrypt(nonce, ciphertext, b"")

    msg_timestamp, message = re.split(b"__", message, maxsplit=1)
    msg_timestamp = int.from_bytes(msg_timestamp, byteorder="big")
    timestamp = calendar.timegm(dt.now().utctimetuple())
    if (timestamp - msg_timestamp) > timeout:
        raise DecryptionTimeout("The message has an expired timeout")
    return message.decode("utf-8")


def handle_decrypt_message(args):
    key = base64.b64decode(args.key.encode("utf-8"))
    ciphertext = base64.b64decode(args.message.encode("utf-8"))
    message = decrypt_message(key, ciphertext, timeout=args.timeout)
    print(str(message), file=args.outfile)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outfile", default=sys.stdout, type=argparse.FileType("w"), help="The output file")

    subparsers = parser.add_subparsers(help="Available sub-commands")

    gen_parser = subparsers.add_parser("generate_key", help="Generate a random key string")
    gen_parser.set_defaults(func=handle_generate_key)

    encrypt_parser = subparsers.add_parser("encrypt", help="Encrypt a message")
    encrypt_parser.add_argument("--key", required=True, help="The encryption key")
    encrypt_parser.add_argument("--message", required=True, help="Message to encrypt")
    encrypt_parser.set_defaults(func=handle_encrypt_message)

    decrypt_parser = subparsers.add_parser("decrypt", help="Decyrpt a message")
    decrypt_parser.add_argument("--key", required=True, help="The encryption key")
    decrypt_parser.add_argument("--message", required=True, help="Message to decrypt")
    decrypt_parser.add_argument(
        "--timeout",
        default=MESSAGE_TIMEOUT,
        type=int,
        help="A message timeout. Decryption will fail for older messages",
    )
    decrypt_parser.set_defaults(func=handle_decrypt_message)

    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    args.func(args)
