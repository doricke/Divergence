
require 'input_file.rb'

class TaxonomyDiverge

attr_accessor :q_ave
attr_accessor :q_max

################################################################################
# Copyright (C) 2015  Darrell O. Ricke, PhD
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
################################################################################

################################################################################
def initialize
  @q_ave = {}
  @q_max = {}
end  # initialize

################################################################################
def read_divergence( taxfile )
  # puts "read_divergence called: #{taxfile}"
  inf = InputFile.new( taxfile )
  inf.open_file()
  # line = inf.next_line()   # header line
  while ( inf.is_end_of_file? == false )
    line = inf.next_line()
    if ( inf.is_end_of_file? == false )
      tokens = line.split( "," )
      name = tokens[2]
      q_ave_div = tokens[5]
      q_max_div = tokens[6]
      n = tokens[10].to_i
      puts "name: #{name}, q1 #{q_ave_div}, q2 #{q_max_div}, n: #{n}"
      if ( n > 2 )
        @q_ave[ name ] = [] if @q_ave[ name ].nil?
        @q_max[ name ] = [] if @q_max[ name ].nil?
        q1 = q2 = 0.0
        q1 = q_ave_div.to_f if (q_ave_div != "NaN")
        a2 = q_max_div.to_f if (q_max_div != "NaN")
        @q_ave[ name ] << q1 if (q1 > 0.0)
        @q_max[ name ] << q2 if (q2 > 0.0)
      end  # if
    end  # if
  end  # while
  inf.close_file()
end  # read_divergence

################################################################################
def snap
  # puts "snap called"
  @q_ave.each do |name, values|
    i = (values.size + 1) / 2
    j = (@q_max[name].size + 1) / 2
    puts "#{name}\t#{(values.sort)[i]}\t#{(@q_max[name].sort)[j]}\t#{values.size}\t#{values}" if values.size > 3
  end  # do
end  # snap

################################################################################

end  # class Msa


################################################################################
def test_divergence( tax_file )
  puts "usage: ruby taxonomy_divergence.rb Vert_tax.csv\n" if ARGV.length < 1
  seqs_file = ARGV[0] if ARGV.length >= 1

  app = TaxonomyDiverge.new
  app.read_divergence( tax_file )
  app.snap
end  # method test_divergence


################################################################################
  test_divergence( "Vert_tax.csv" )
