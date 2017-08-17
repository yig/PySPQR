'''
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
Description: An `asarray` function that wraps a cffi pointer in a numpy.array.
URL: https://gist.github.com/yig/77667e676163bbfc6c44af02657618a6
'''

from __future__ import print_function, division, absolute_import

import numpy

## Create the dictionary mapping ctypes to numpy dtypes.
ctype2dtype = {}

## Integer types
for prefix in ( 'int', 'uint' ):
    for log_bytes in range( 4 ):
        ctype = '%s%d_t' % ( prefix, 8*(2**log_bytes) )
        dtype = '%s%d' % ( prefix[0], 2**log_bytes )
        # print( ctype )
        # print( dtype )
        ctype2dtype[ ctype ] = numpy.dtype( dtype )

## Floating point types
ctype2dtype[ 'float' ] = numpy.dtype( 'f4' )
ctype2dtype[ 'double' ] = numpy.dtype( 'f8' )

# print( ctype2dtype )

def asarray( ffi, ptr, length ):
    ## Get the canonical C type of the elements of ptr as a string.
    T = ffi.getctype( ffi.typeof( ptr ).item )
    # print( T )
    # print( ffi.sizeof( T ) )

    if T not in ctype2dtype:
        raise RuntimeError( "Cannot create an array for element type: %s" % T )

    return numpy.frombuffer( ffi.buffer( ptr, length*ffi.sizeof( T ) ), ctype2dtype[T] )

def test():
    from cffi import FFI
    ffi = FFI()

    N = 10
    ptr = ffi.new( "float[]", N )

    arr = asarray( ffi, ptr, N )
    arr[:] = numpy.arange( N )

    for i in range( N ):
        print( arr[i], ptr[i] )

if __name__ == '__main__':
    test()
