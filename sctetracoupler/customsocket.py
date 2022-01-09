import struct
import socket
import select
import time
import subprocess
import numpy as np

class CustomSocketClient(object):
    """Low level details of the socket communication"""
    def __init__(self,tcp_id='127.0.0.1', tcp_port=3333):
        self.tcp_id = tcp_id
        self.tcp_port = tcp_port
        self.handshake = 178278912
        
    def start(self):
        """() -> None. Open the low level socket connection. Blocks but allows the Python thread scheduler to run.

        """
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect((self.tcp_id,self.tcp_port))
        self.conn = self.socket
        self.send_data(self.handshake)

    def send_data(self, value):
        """(value: any) -> None. Send value to Itasca software. value must be int, float, length two list of doubles, length three list of doubles or a string.

        """
        while True:
            _, write_ready, _ = select.select([], [self.conn], [], 0.0)
            if write_ready: break
            else: time.sleep(1e-8)

        if type(value) == int:
            self.conn.sendall(struct.pack("i", 1))
            self.conn.sendall(struct.pack("i", value))
        elif type(value) == float:
            self.conn.sendall(struct.pack("i", 2))
            self.conn.sendall(struct.pack("d", value))
        elif type(value) == list and len(value)==2:
            float_list = [float(x) for x in value]
            self.conn.sendall(struct.pack("i", 5))
            self.conn.sendall(struct.pack("dd", float_list[0], float_list[1]))
        elif type(value) == list and len(value)==3:
            float_list = [float(x) for x in value]
            self.conn.sendall(struct.pack("i", 6))
            self.conn.sendall(struct.pack("ddd", float_list[0],
                                          float_list[1], float_list[2]))
        elif type(value) == str:
            length = len(value)
            self.conn.sendall(struct.pack("ii", 3, length))
            buffer_length = 4*(1+(length-1)/4)
            format_string = "%is" % buffer_length
            value += " "*(buffer_length - length)
            self.conn.sendall(struct.pack(format_string, value))
        else:
            raise Exception("unknown type in send_data")

    def wait_for_data(self):
        """() -> None. Block until data is available. This call allows the Python thread scheduler to run.

        """
        while True:
            input_ready, _, _ = select.select([self.conn],[],[], 0.0)
            if input_ready: return
            else: time.sleep(1e-8)

    def read_type(self, type_string):
        """(type: str) -> any. This method should not be called directly. Use the read_data method.

        """
        byte_count = struct.calcsize(type_string)
        bytes_read = 0
        data = ''
        self.wait_for_data()
        while bytes_read < byte_count:
            self.wait_for_data()
            data_in = self.conn.recv(byte_count - bytes_read)
            data += data_in
            bytes_read += len(data)
        assert len(data)==byte_count, "bad packet data"
        return data

    def read_data(self):
        """() -> any. Read the next item from the socket connection."""
        raw_data = self.read_type("i")
        type_code, = struct.unpack("i", raw_data)
        if type_code == 1:     # int
            raw_data = self.read_type("i")
            value, = struct.unpack("i", raw_data)
            return value
        elif type_code == 2:   # float
            raw_data = self.read_type("d")
            value, = struct.unpack("d", raw_data)
            return value
        elif type_code == 3:   # string
            length_data = self.read_type("i")
            length, = struct.unpack("i", length_data)
            buffer_length = (4*(1+(length-1)/4))
            format_string = "%is" % buffer_length
            data = self.read_type(format_string)
            return data [:length]
        elif type_code == 5:   # V2
            raw_data = self.read_type("dd")
            value0, value1 = struct.unpack("dd", raw_data)
            return [value0, value1]
        elif type_code == 6:   # V3
            raw_data = self.read_type("ddd")
            value0, value1, value3 = struct.unpack("ddd", raw_data)
            return [value0, value1, value3]
        assert False, "Data read type error %i" % type_code

    def close(self):
        """() -> None. Close the active socket connection."""
        self.conn.close()