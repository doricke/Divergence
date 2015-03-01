
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

class TaxonomyLevel

###############################################################################

# Average percent identities
attr_accessor :averages
# Taxonomy level
attr_accessor :level
# Maximum mismatches observed
attr_accessor :max_mismatches
# Maximum number of aligned residues observed
attr_accessor :max_residues
# Names
attr_accessor :names
# Common taxonomy
attr_accessor :taxonomy
# Number of sequences/species included in the current MSA.
attr_accessor :total
# Total number of residue mismatches
attr_accessor :total_mismatches
# Total number of residues
attr_accessor :total_residues


###############################################################################
def initialize
  @averages = {}
  @level = ""
  @max_mismatches = 0
  @max_residues = 0
  @names = {}
  @taxonomy = ""
  @total_mismatches = 0  
  @total_residues = 0
end  # initialize

###############################################################################
def add_comparison( identities, mismatches, residues, name1, name2 )
  @max_mismatches = mismatches if (mismatches > @max_mismatches)
  @max_residues = residues if (residues > @max_residues)
  @total_mismatches += mismatches
  @total_residues += residues
  @names[ name1 ] = true
  @names[ name2 ] = true

  # Track the percent identity averages.
  return if residues < 1
  percent = (identities.to_f * 100.0) / residues.to_f
  @averages[ name1 ] = [] if @averages[ name1 ].nil?
  @averages[ name2 ] = [] if @averages[ name2 ].nil?
  @averages[ name1 ] << percent
  @averages[ name2 ] << percent
end  # add_comparison

###############################################################################
def calculate_v
  return (1.0 - fraction) / 0.451683061
end  # calculate_v

###############################################################################
def flag_names
  names = []
  total = 0.0
  count = 0
  @averages.each do |name, aves|
    aves.each do |ave|
      total += ave
      count += 1
    end  # do
  end  # do
  average = total / count.to_f

  # Flag outliers.
  @averages.each do |name, aves|
    total = 0.0
    aves.each do |ave|
      total += ave
    end  # do
    ave = total / aves.size.to_f

    # flag outliers.
    names << name if ( ave + 7.0 < average ) 
  end  # do

  return names.join( "," )
end  # flag_names

###############################################################################
def fraction
  return (@total_residues.to_f - @total_mismatches.to_f) / @total_residues.to_f
end  # fraction

###############################################################################
def get_names
  return @names.keys
end  # get_names

###############################################################################
def mismatches
  # return @max_mismatches
  total = @names.keys.size
  return 0 if (total < 1)
  return @total_mismatches / total 
end  # mismatches

###############################################################################
def percent 
  return 0 if (@names.keys.size < 1) || (@total_residues < 1)
  return ((@total_residues - @total_mismatches) * 100) / @total_residues
end  # percent

###############################################################################
def residues
  total = @names.keys.size
  return 0 if (total < 1)
  return @total_residues / total
end  # residues

###############################################################################

end  # TaxonomyLevel

