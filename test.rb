require 'nmatrix'
require 'ruby-prof'
require 'benchmark'

num = 100000
n = NMatrix.new([num], [1]*num)
a = [1]*num

#result = RubyProf.profile do 
#  n.sum
#end
#
#printer = RubyProf::GraphHtmlPrinter.new(result)
#printer.print(STDOUT)
#
#exit

Benchmark.bm do |x|
    x.report("nm") { n.sum }
    x.report("arry") {a.inject(:+)}
end
