require 'nmatrix'
require 'pp'

dtype = :float32
a = NMatrix.new(3, [-2,4,-3,3,-2,1,0,-4,3], dtype: dtype)

#this is the a and ipiv that we should get from getrf, but getrf doesn't work currently 
a_intermediate = NMatrix.new([3,3],[-3.00, 0.00, -10.00,-1.0/3,1.00,8.0/3,-1.0/3,0.0,2.0/3], dtype: dtype)
ipiv = [3,2,3]

NMatrix::LAPACK::clapack_getri(:row, 3, a_intermediate, 3, ipiv)

puts "calculated result:"
pp a_intermediate
puts "expected result:"
puts "[
  [-5.0, -0.0, -2.0]   [-4.0,  1.0, -1.0]   [ 1.5,  0.0,  0.5] ]"

