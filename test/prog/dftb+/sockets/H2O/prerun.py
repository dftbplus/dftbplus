#!/usr/bin/env python
from __future__ import print_function
import struct
import socket
import os
import numpy as np
import numpy.linalg as la
from sockettools import frac2cart, readgen, receive_all, a0

# Expecting two geometry steps of communication with DFTB+
NR_STEPS = 2

def connect():
    server_address = '/tmp/ipi_dftb'

    # Make sure the socket does not already exist
    try:
        os.unlink(server_address)
    except OSError:
        if os.path.exists(server_address):
            raise

    serversocket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    serversocket.bind(server_address)
    serversocket.listen(1) # become a server socket, maximum 1 connections
    connection, address = serversocket.accept()

    return connection



def main():

    specienames, species, coords, origin, latvecs = readgen("ice.gen")
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

    connection.shutdown(socket.SHUT_RDWR)
    connection.close()


if __name__ == '__main__':
    main()
