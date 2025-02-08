#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

from __future__ import print_function
import struct
import socket
import os
import numpy as np
import numpy.linalg as la
from sockettools import frac2cart, readgen, receive_all, a0

# Expecting five geometry steps of communication with DFTB+, sending
# the same structure each time
NR_STEPS = 5

def connect():

    pid = os.getpid()

    server_address = '/tmp/ipi_dftb%i' % pid

    # Make sure the socket does not already exist
    try:
        os.unlink(server_address)
    except OSError:
        if os.path.exists(server_address):
            raise

    # write file for dftb_in.hsd to include:
    file = open("file.hsd","w")
    file.write('# The externally set filename for this run\n')
    file.write("+Driver = +Socket {\n")
    file.write('  !File = "dftb%i"\n' % pid)
    file.write("}\n")
    file.close()
    # plain text file with the same information
    file = open("file.txt","w")
    file.write("/tmp/ipi_dftb%i" % pid)
    file.close()

    serversocket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    serversocket.bind(server_address)
    serversocket.listen(1) # become a server socket, maximum 1 connections
    connection, address = serversocket.accept()

    return connection



def main():

    specienames, species, coords, origin, latvecs = readgen("geo.gen")
    coords /= a0
    if (latvecs):
        latvecs /= a0
        tPeriodic = 1
    else:
        latvecs = np.empty((3, 3), dtype=float)
        tPeriodic = None

    connection = connect()

    for iStep in range(NR_STEPS+1):
        print("Step %i" % iStep)
        connection.sendall(b'POSDATA     ')
        connection.sendall(latvecs)

        if tPeriodic:
           connection.sendall(la.inv(latvecs))
        else:
           connection.sendall(np.empty((3, 3), dtype=float))

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
        if tPeriodic:
            print(unpacked_data)
            # contains numerical junk if not peridic, so do nothing
            # otherwise

        # dummy '0' at the end of the data
        buf = receive_all(connection, 4)
        unpacked_data = struct.unpack('i', buf)
        print(unpacked_data)

    connection.shutdown(socket.SHUT_RDWR)
    connection.close()


if __name__ == '__main__':
    main()
