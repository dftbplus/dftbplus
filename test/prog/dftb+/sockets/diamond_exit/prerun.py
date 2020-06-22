#!/usr/bin/env python3
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
    serversocket.bind(('localhost', 0))
    serversocket.listen(1)
    port = serversocket.getsockname()[1]
    # write file for dftb_in.hsd to include:
    file = open("port.hsd","w")
    file.write('# The externally set port number for this run\n')
    file.write("+Driver = +Socket {\n")
    file.write("  !Port = %i\n" % port)
    file.write("}\n")
    file.close()
    # plain text file with the same information
    file = open("port.txt","w")
    file.write("%i" % port)
    file.close()
    connection, address = serversocket.accept()
    return connection


def main():

    specienames, species, coords, origin, latvecs = readgen("diamond.gen")
    coords /= a0
    latvecs /= a0

    connection = connect()

    for iStep in range(NR_STEPS):
        print("Step %i" % iStep)
        connection.sendall(b'POSDATA     ')
        connection.sendall(latvecs)

        connection.sendall(la.inv(latvecs))
        connection.sendall(np.int32(len(coords)))
        connection.sendall(coords)

        connection.sendall(b'GETFORCE    ')
        # needs work:
        buf = receive_all(connection, 12)
        if (buf != b'FORCEREADY  '):
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
    connection.sendall(b'EXIT        ')

    time.sleep(2)
    connection.shutdown(socket.SHUT_RDWR)
    connection.close()


if __name__ == '__main__':
    main()
