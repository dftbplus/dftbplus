#!/usr/bin/env python
from __future__ import print_function
import struct
import socket
import time
import numpy as np
import numpy.linalg as la
from sockettools import frac2cart, readgen, receive_all, a0

# Expecting only one geometry step of communication with DFTB+
NR_STEPS = 1

def connect():
    serversocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    serversocket.bind(('localhost', 21013))
    serversocket.listen(1)
    connection, address = serversocket.accept()
    return connection


def main():

    specienames, species, coords, origin, latvecs = readgen("diamond.gen")
    coords /= a0
    latvecs /= a0

    connection = connect()

    for iStep in range(NR_STEPS):
        print("Step %i" % iStep)
        connection.sendall('POSDATA     ')
        connection.sendall(latvecs)

        connection.sendall(la.inv(latvecs))
        connection.sendall(np.int32(len(coords)))
        connection.sendall(coords)

        connection.sendall('GETFORCE    ')
        # needs work:
        buf = receive_all(connection, 12)
        if (buf != 'FORCEREADY  '):
            raise ValueError('Unexpected value of "GETFORCE": "%s"!' % buf)

        # expecting energy and number of atoms
        buf = receive_all(connection, 12)
        unpacked_data = struct.unpack('di', buf)
        print(unpacked_data)

        # forces
        buf = receive_all(connection, len(coords)*3*8)
        frmat = '%i' % (3*len(coords))
        frmat += 'd'
        unpacked_data = struct.unpack(frmat, buf)
        print(unpacked_data)

        # stress x lattive volume
        buf = receive_all(connection, 9*8)
        unpacked_data = struct.unpack('9d', buf)
        print(unpacked_data)

        # dummy '0' at the end of the data
        buf = receive_all(connection, 4)
        unpacked_data = struct.unpack('i', buf)
        print(unpacked_data)


    print("Requesting clean shutdown of DFTB+")
    connection.sendall('EXIT        ')

    time.sleep(2)
    connection.shutdown(socket.SHUT_RDWR)
    connection.close()


if __name__ == '__main__':
    main()
