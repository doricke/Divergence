
require 'fasta_iterator.rb'
require 'fasta_sequence.rb'
require 'tuples.rb'
require 'output_file.rb'
require 'taxonomy_level.rb'

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
# This class represents a protein multiple sequence alignment.
class Diverge
         
################################################################################

# Minimum percent identity to consider before sequence is considered too diverged
MIN_PERCENT_IDENTITY = 40

# Limit for number of adjacent spacers before excluding from calculations.
SPACER_LIMIT = 15
         
################################################################################

# Final MSA alignment
attr_accessor :aligned

# Bird residues in MSA
attr_accessor :aves

# human gene chromosome position
attr_accessor :chromosome

# Number of excluded spacer residues
attr_accessor :excluded_spaces

# Number of invariant residues
attr_accessor :invariant

# CpG residues
attr_accessor :is_cpg

# Mammal residues in MSA
attr_accessor :eutheria

# FASTA file name
attr_accessor :fasta_name

# Longest protein sequence
attr_accessor :longest

# gene chromosome map position for human gene
attr_accessor :map

# Protein multiple sequence alignment hash - hashed by sequence name
attr_accessor :msa

# Observed residues at each position.
attr_accessor :observed

# Estimated q value
attr_accessor :q

# Taxonomy leaves
attr_accessor :taxonomy_leaves

# Taxonomy tree
attr_accessor :taxonomy_tree

# Word tuples
attr_accessor :tuples

# Sequences
attr_accessor :seqs

# Number of adjacent spacer residues
attr_accessor :spacers

# Number of variable positions in the current MSA (positions with gaps in first sequence not counted).
attr_accessor :variable


################################################################################
def initialize
  @aligned = {}
  @aves = {}
  @chromosome = {}
  @excluded_spacers = 0
  @invariant = 0
  @is_cpg = {}
  @eutheria = {}
  @longest = 0
  @map = {}
  @msa = {}
  @observed = []
  @q = 0.9980
  @tuples = {}
  @seqs = {}
  @spacers = []
  @taxonomy_leaves = {}
  @taxonomy_tree = {}
  @variable = 0
end  # initialize


################################################################################
# This method counts up the number of estimated CpG positions in the alignment.
def count_cpg
  cpg_total = 0
  @is_cpg.each do |i, flag|
    cpg_total += 1 if flag
  end  # do
  return cpg_total
end  # count_cpg


################################################################################
# This method counts up the number of excluded spacers in the alignment.
def count_excluded
  spacers_excluded = 0
  names = @seqs.keys()
  first = names[0]
  @spacers.each do |i, count|
    spacers_excluded +=1 if (! count.nil?) && (count >= SPACER_LIMIT) && (@aligned[ first ][ i ] != '.')
  end  # do
  return spacers_excluded
end  # count_excluded


################################################################################
# Tally class-specific residues.
def count_residues( seqs_file )
  out = OutputFile.new( seqs_file + "_class.txt" )
  out.open_file
  count = 1
  @aligned.each do |name, segment|
    tax = @seqs[name].annotation["taxonomy"]
    organism = @seqs[ name ].common_name
    if tax.include?( "Mammalia" )
      out.write( "#{count}\t#{name}\t#{organism}\tMammal\n" )
      count += 1
    end  # if
  end  # do

  @aligned.each do |name, segment|
    tax = @seqs[name].annotation["taxonomy"]
    organism = @seqs[ name ].common_name
    if tax.include?( "Aves" )
      out.write( "#{count}\t#{name}\t#{organism}\tBird\n" )
      count += 1
    end  # if
  end  # do

  names = @seqs.keys()
  first = names[0]
  pos = 1
  if ! @aligned[ first ].nil?
    for i in 0...@aligned[ first ].size do
      # Check for any residues.
      residues = 0
      @aligned.each do |name, segment|
        tax = @seqs[name].annotation["taxonomy"]
        residues += 1 if (tax.include?( "Mammalia" ) || tax.include?( "Aves" ) ) && (@aligned[ name ][ i ] != '.')
      end  # do
  
      if ( residues > 0 )
        out.write( "#{pos} " )
        @aligned.each do |name, segment|
          tax = @seqs[name].annotation["taxonomy"]
          if tax.include?( "Mammalia" )
            out.write( @aligned[ name ][ i ] )
          end  # if
        end  # do
    
        out.write( " " )
        @aligned.each do |name, segment|
          tax = @seqs[name].annotation["taxonomy"]
          if tax.include?( "Aves" )
            out.write( @aligned[ name ][ i ] )
          end  # if
        end  # do
       
        out.write( "\n" )
        pos += 1
      end  # if
    end  # for
  end  # if
  
  out.close_file
end  # count_residues


################################################################################
# Count the variable residues in the MSA.
def count_variable
  @variable = 0
  names = @seqs.keys()
  first = names[0]
  return if @seqs[ first ].nil?

  # Traverse each sequence.
  @aligned.each do |name, segment|
    tax = @seqs[name].annotation["taxonomy"]
    # if tax.include?( "Euteleostomi" )
    if tax.include?( "Mammalia" )
      # puts "Considering: #{name}\t#{tax}"

      # Traverse each position in the alignment segment.
      for i in 0...segment.size do
        @observed[i] = {} if @observed[i].nil?
        if (!@aligned[first][i].nil? ) && (@aligned[first][i] != '.')
          residue = @aligned[name][i]
          if ( residue != '.' ) && ( residue != 'X' ) && ( residue != '-' )
            @observed[i][residue] = 0 if @observed[i][residue].nil?
            @observed[i][residue] += 1
          end  # if
        end  # if
      end  # for
    end  # if

    if tax.include?( "Aves" )
      # Traverse each position in the alignment segment.
      for i in 0...segment.size do
        @aves[i] = {} if @aves[i].nil?
        if (!@aligned[first][i].nil? ) && (@aligned[first][i] != '.')
          residue = @aligned[name][i]
          if ( residue != '.' ) && ( residue != 'X' ) && ( residue != '-' )
            @aves[i][residue] = 0 if @aves[i][residue].nil?
            @aves[i][residue] += 1
          end  # if
        end  # if
      end  # for
    end  # if      
  end  # do

  for i in 0...@observed.size do
    if (@observed[i].keys.size > 1)
      @variable += 1 
    else
      @invariant += 1 if (@observed[i].keys.size == 1 )
    end  # if
    # puts "#{i}\t#{@observed[i]}"
  end  # for
 
  puts "Variable: #{@variable}"
end  # count_variable
  
  
################################################################################
def is_cpg_position( counts )
  return false if counts.size <= 1
      
  # puts "#{counts} counts.size: #{counts.size}"
      
  total = 0
  residues = {}
  residues["R"] = residues["P"] = residues["S"] = residues["T"] = residues["V"] = 0
  residues["A"] = residues["D"] = residues["Q"] = residues["G"] = 0
      
  counts.each do |residue, count|
    case residue
      when "R", "r", "H", "h", "Q", "q", "C", "c", "W", "w", "*"
        residues["R"] += count
      when "P", "p", "L", "l"
        residues["P"] += count
      # when "S", "s", "L", "l"
        #  residues["S"] += count
      # when "T", "t", "M", "m"
        #  residues["T"] += count
      # when "V", "v", "I", "i", "M", "m"
        #   residues["V"] += count
      # when "A", "a", "T", "t"
        #   residues["A"] += count
      # when "D", "d", "N", "n"
        #   residues["D"] += count
      # when "Q", "q", "K", "k"
        #   residues["Q"] += count
      # when "G", "g", "S", "s", "R", "r"
        #   residues["G"] += count
    end  # case
    total += count
  end  # do
      
  # Check for CpG position
  residues.each do |residue, count|
    if ( (count == total) && (! counts[residue].nil?) )
      # puts "--: residue: #{residue} count:#{count} total #{total}"
      return true
    end  # if
  end  # do
  return false
end  # is_cpg_position
  
################################################################################
def determine_cpg
  for i in 0...@observed.size do
    if (@observed[i].keys.size > 1)
      if is_cpg_position( @observed[i] )
        @is_cpg[i] = true
        puts "CpG: #{i}\t#{@observed[i]}"
      end
    else
      @is_cpg[i] = false
    end  # if
  end  # for
end  # determine_cpg
      
      
################################################################################
def determine_spacers
          
  # Count up the number of adjacent spacer residues
  j = 0
  while ( j < @observed.size) do
    count = 1
    if (@observed[j].keys.size > 1)
      i = j + 1
      while ( ( i < @observed.size) && (@observed[i].keys.size > 1) ) do
        i += 1
        count += 1
      end  # while
                  
      for k in j...i do
        @spacers[k] = count
      end  # do
                  
      if (i-j+1 >= 10)
        puts "Spacers: #{j+1}-#{i} count #{(i-j+1)}" 
        @excluded_spacers += (i-j+1)
      end  # if
    end  # if
    j += count
  end  # while
                  
end  # determine_spacers


################################################################################
def show_adaptation
  for i in 0...@observed.size do
    if (! @aves[i].nil?)
      if ( (@observed[i].size == 1) || (@aves[i].size == 1) )  &&
         # ( @observed[i].keys.sort != @aves[i].keys.sort )
         ( ( @observed[i].keys & @aves[i].keys ).empty? )
        puts "#{@observed[i]} birds: #{@aves[i]}" if (@aves[i].size >= 5)
      end  # if
    end  # if
  end  # for
end  # show_adaptation


################################################################################
# This method determines the consensus for a position.
def consensus
  cons = ""
  return cons
end  # method consensus


################################################################################
# This method reads in the FASTA sequences from the file.
def load_fastas( fasta_filename )
  @fasta_name = fasta_filename
  in_fasta = FastaIterator.new( fasta_filename )
  in_fasta.open_file
  @seqs = in_fasta.read_fastas
  @seqs.each do |name, seq|
    seq.parse_annotation
    len = seq.sequence_data.size
    @longest = len if (len > @longest)
  end  # do
  return @seqs
end  # method load_fastas


################################################################################
def calculate_tuples
  @seqs.each do |name, fasta|
    tuples = Tuples.new
    tuples.calculate_words( fasta.sequence_data )
    @tuples[ fasta.sequence_name ] = tuples
  end  # do
end  # calculate_tuples


################################################################################
def compare_tuples
  names = @seqs.keys()
  first = names[0]
  return if @seqs[ first ].nil?
  seq1 = @seqs[ first ].sequence_data
  @msa[ first ] = @seqs[ first ].sequence_data
  for name in names[1..-1] do
    # puts "Aligning: #{name}"
    seq2 = @seqs[ name ].sequence_data
    @tuples[ first ].calculate_common( @tuples[ name ].words, @tuples[ name ].words_positions )
    # @tuples[ first ].show_common_words( @tuples[ name ].words_positions )
    @msa[ name ] = @tuples[ first ].align_words( seq1, seq2, @tuples[ name ].words_positions, name )
    # puts
  end  # for
end  # compare_tuples


################################################################################
# This method drops in the insert sequence and gap fills to the width for this position.
def add_extra( ins_seq, width )
  seq = ins_seq
  count = ins_seq.length
  while ( count < width ) do
    seq << "."
    count += 1
  end  # while
  return seq
end  # add_extra


################################################################################
# This method adds the insert for this sequence or gaps where needed.
def add_inserts( align, inserts, widths, name )
  return nil if align.nil? || ( align.keys.size < 1 )
  seq = ""
  for i in 1..align.keys.sort.last do
    if (! align[i].nil?)
      seq << align[i]
    else
      seq << "."
    end  # if

    ins_seq = ""
    ins_seq = inserts[i][name] if (! inserts[i].nil?) && (! inserts[i][name].nil?)
    seq << add_extra( ins_seq, widths[i] ) if (! widths[i].nil?)
  end  # for

  return seq
end  # add_inserts


################################################################################
def compare_aligned( align1, align2 )
  return 0, 0, 0 if align1.nil? || align2.nil? || (align2.size < 1)
  identities = 0
  mismatches = 0
  compared = 0
  for i in 0...align1.size do
    @spacers[i] = 0 if @spacers[i].nil?
    @is_cpg[i] = false if @is_cpg[i].nil?
    if (i < align2.size) && ( @spacers[i] < SPACER_LIMIT ) && ( @is_cpg[i] == false )
      if ( align1[i] != '.' ) && ( align2[i] != '.' ) &&
         ( align1[i] != 'X' ) && ( align2[i] != 'X' ) &&
         ( align1[i] != '-' ) && ( align2[i] != '-' )
        compared += 1
        if ( align1[i] == align2[i] )
          identities += 1
        else
          mismatches += 1
        end  # if
      end  # if
    end  # if
  end  # for
  return identities, mismatches, compared
end  # compare_aligned


################################################################################
# This method determines the insert maximum size at each position.
def size_inserts( inserts )
  widths = {}
  inserts.each do |pos, extras| 
    extras.each do |name, ins_seq|
      ins_size = ins_seq.size
      if ( widths[ pos ].nil? )
        widths[ pos ] = ins_size
      else
        # Check if next insert is larger than the previous at this position.
        widths[ pos ] = ins_size if ( ins_size > widths[ pos ] )
      end
    end  # do
  end  # do

  return widths
end  # size_inserts


################################################################################
def calculate_q( mismatches, v, residues, n )
  if ( v > 0.0 ) && ( mismatches < v * residues ) && ( n > 0 )
    q = Math.exp( Math.log( 1.0 - mismatches /(v * residues)) / n )
    q = 0.998 if q.nil?
    return q
  end  # if

  return 0.998  # default q
end  # calculate_q


################################################################################
def calculate_n( mismatches, v, residues, q )
  if ( v > 0.0 ) && ( mismatches < v * residues ) && (! q.nil?)
    n = Math.log(1.0 - mismatches/(v * residues )) / Math.log( q )
    return n
  end  # if

  return 0.0
end  # calculate_n


################################################################################
def calculate_n_range( taxonomy_level, v, q_ave, q_max )
  residues = taxonomy_level.residues
  mismatches = taxonomy_level.mismatches
  max_residues = taxonomy_level.max_residues
  max_mismatches = taxonomy_level.max_mismatches
  n_ave = calculate_n( mismatches, v, residues, q_ave )

  max_residues = taxonomy_level.max_residues
  max_mismatches = taxonomy_level.max_mismatches
  n_max = calculate_n( max_mismatches, v, max_residues, q_max )

  return n_ave, n_max
end  # calculate_n_range


################################################################################
def calculate_q_range( taxonomy_level, v, n_ave, n_max )
  residues = taxonomy_level.residues
  mismatches = taxonomy_level.mismatches
  max_residues = taxonomy_level.max_residues
  max_mismatches = taxonomy_level.max_mismatches
  q_ave = calculate_q( mismatches, v, residues, n_ave )

  max_residues = taxonomy_level.max_residues
  max_mismatches = taxonomy_level.max_mismatches
  q_max = calculate_q( max_mismatches, v, max_residues, n_max )

  return q_ave, q_max
end  # calculate_q_range


################################################################################
def add_taxonomy_leaves( level, n_ave, n_max, names, depth )
  species = names.split( "|" )
  species.each do |organism|
    @taxonomy_leaves[ organism ] = {} if @taxonomy_leaves[ organism ].nil?
    if (@taxonomy_leaves[ organism ][ "node" ].nil? ) || (depth > @taxonomy_leaves[ organism ][ "depth" ])
      @taxonomy_leaves[ organism ][ "node" ] = level
      @taxonomy_leaves[ organism ][ "n_ave" ] = n_ave
      @taxonomy_leaves[ organism ][ "n_max" ] = n_max
      @taxonomy_leaves[ organism ][ "name" ] = organism
      @taxonomy_leaves[ organism ][ "depth" ] = depth
    end  # if
  end  # do
end  # add_taxonomy_leaves


################################################################################
def add_taxonomy_node( taxonomy_tree, taxonomy, level, n_ave, n_max, names, depth )
  tax = taxonomy.split( "|" )
  top_level = tax[ 0 ]
  taxonomy_tree[ top_level ] = {} if ( taxonomy_tree[ top_level ].nil? )
  if ( level == top_level )
    taxonomy_tree[ level ][ "name" ] = level
    taxonomy_tree[ level ][ "leaves" ] = names
    taxonomy_tree[ level ][ "n_ave" ] = n_ave
    taxonomy_tree[ level ][ "n_max" ] = n_max
    add_taxonomy_leaves( level, n_ave, n_max, names, depth )
  else
    add_taxonomy_node( taxonomy_tree[ top_level ], tax[ 1..-1 ].join( "|" ), level, n_ave, n_max, names, depth+1 )
  end  # if
end  # add_taxonomy_node


################################################################################
def write_spaces( out, spaces )
  for i in 0...spaces do
    out.write( " " )
  end  # for
end  # write_spaces


################################################################################
def write_children( out, spaces, level )
  write_spaces( out, spaces )
  out.write( "\"children\": [\n" )
  count = 0
  @taxonomy_leaves.each do |name, values|
    # puts "level: #{level}, node: #{values['node']}"
    if ( values[ "node" ] == level )
      out.write( "," ) if (count > 0)
      write_spaces( out, spaces+1 )
      size = (values[ "n_ave" ] + values[ "n_max" ]) / 2.0
      out.write( "{\"name\": \"#{name}\", \"size\": #{size/2.0}}\n" )
      count += 1
    end  # if
  end  # do
end  # write_children


################################################################################
def write_json_nodes( out, spaces, taxon_tree, previous )
  return if taxon_tree.keys.size < 1

  taxon_tree.each do |name, node|
    if taxon_tree[ name ].is_a?(::Hash) && (! taxon_tree[ name ][ "name" ].nil?)
      out.write( "," ) if previous
      write_spaces( out, spaces )
      out.write( "{\n" )
      write_spaces( out, spaces )
      out.write( " \"name\": \"#{name}\",\n" )
      # Scan for children
      write_children( out, spaces, name )
      write_json_nodes( out, spaces+1, node, false ) if node.is_a?(::Hash)
      write_spaces( out, spaces+1 )
      out.write( "]\n" )
      write_spaces( out, spaces )
      out.write( "}\n" )
    else
      write_json_nodes( out, spaces, node, false ) if node.is_a?(::Hash)
    end  # if

    previous = true if node.is_a?(::Hash)
  end  # do
end  # write_json_nodes

################################################################################
def write_taxonomy_json( filename )
  out = OutputFile.new( filename + ".json" )
  out.open_file
  spaces = 0
  previous = false
  write_json_nodes( out, spaces, @taxonomy_tree, previous )
  out.close_file
end  # write_taxonomy_json

################################################################################
def get_leaves( level )
  children = ""
  @taxonomy_leaves.each do |organism, leaf|
    if ( leaf[ "node" ] == level )
      children += "," if (children.size > 0)
      children += "#{leaf['name']}"
      children += ":#{(leaf['n_ave']/2.0+0.5).to_i}" if (! leaf['n_ave'].nan?) && (! leaf['n_ave'].infinite?)
    end  # if
  end  # do 

  puts "*** Node: #{level} terminal leaves: #{children}"
  return children
end  # get_leaves

################################################################################
def newick_nodes( taxon_tree, node_name )
  puts "newick_nodes: #{node_name} #{taxon_tree}"

  return "" if taxon_tree.keys.size < 1

  str = ""
  str = "(" if (node_name.size > 0)
  count = 0
  names = taxon_tree.keys
  for i in 0...names.size do
    name = names[i]
    puts "  name: #{name}"
    if taxon_tree[ name ].is_a?(::Hash)
      child_str = ""
      child_str += newick_nodes( taxon_tree[ name ], name ) if ! taxon_tree[ name ].nil?
      if ( child_str.size > 0 )
        str += "," if (count > 0)
        str += child_str
        count += 1
      end  # if
    else
      if ( name != "name" ) && ( name != "leaves" ) && ( name != "n_ave" ) && ( name != "n_max" ) && ( name != "depth" )
        children = get_leaves( taxon_tree[ name ] )
        if ( children.size > 0 )
          str += "," if (count > 0)
          str += children
          count += 1
        # else
        #   str += name
        #   count += 1
        end  # if
      else
        if ( name == "name" )
          children = get_leaves( taxon_tree[ name ] )
          if ( children.size > 0 )
            str += "," if (count > 0)
            str += children
            count += 1
          # else
          #   str += taxon_tree[ "name" ]
          #   str += ":#{(taxon_tree[ "n_ave" ]/2.0+0.5).to_i}" if (! taxon_tree[ "n_ave" ].nil?)
          #   count += 1
          end  # if

        end  # if
      end  # if
    end  # if
  end  # for

  if (node_name.size > 0)
    str += ")#{node_name}" 
    str += ":#{(taxon_tree['n_ave']/2.0+0.5).to_i}" if (! taxon_tree[ 'n_ave' ].nil?) && (! taxon_tree[ 'n_ave' ].nan?) && (! taxon_tree[ 'n_ave' ].infinite?)
  end  # if

  return str
end  # newick_nodes

################################################################################
def write_taxonomy_newick( filename )
  out = OutputFile.new( filename + ".newick" )
  out.open_file
  spaces = 0
  previous = false
  str = newick_nodes( @taxonomy_tree, "" )
  out.write( "#{str}\n" )
  out.close_file
end  # write_taxonomy_newick

################################################################################
def write_comparisons( seqs_file )
  puts "Processing: #{seqs_file}"
  names = @seqs.keys()
  first = names[0]
  return if @seqs[ first ].nil?

  # Set up the sequence taxonomies.
  taxons = {}
  max_length = 0
  @seqs.each do |name, fasta|
    tax = fasta.annotation["taxonomy"]
    taxons[ name ] = []
    taxons[ name ] = tax.split( "; " ) if ! tax.nil?

    # Identify the longest sequence.
    max_length = fasta.sequence_data.size if fasta.sequence_data.size > max_length
  end  # do

  # Don't compute on short sequences.
  return if max_length < 100

  max_mismatches = 0
  residues = 0
  name_i = names[ 0 ]
  alignment_count = 1  # count the number of good alignments
  for j in 1...names.size do
    name_j = names[ j ]
    identities, mismatches, compared = compare_aligned( @aligned[ name_i ], @aligned[ name_j ] )
    percent = 0
    percent = (identities * 100) / compared if compared > 0
    if ( percent < MIN_PERCENT_IDENTITY )
      puts "Reject:\t#{seqs_file}\t#{name_j}" if (compared + 100 >= max_length)
      @aligned[name_j] = ""  # Erase this alignment
    else
      alignment_count += 1 if ( percent > 55 )
    end  # if
    # if (mismatches > max_mismatches) && (compared + 100 >= max_length) && ( percent >= MIN_PERCENT_IDENTITY )
    if (mismatches > max_mismatches) && (compared + 100 >= residues) && ( percent >= MIN_PERCENT_IDENTITY )
      max_mismatches = mismatches 
      residues = compared
    end  # if
  end  # for

  # Return if alignment has less than 6 alignments.
  # return if (alignment_count < 6)

  v = 0.0
  puts "Variable: #{@variable}, residues: #{residues}"
  v = (@variable.to_f / residues.to_f) if (residues > 0)
  v = 0.99 if (v >= 1.0)

  # Compare each of the sequences pairwise.
  out = OutputFile.new( seqs_file + ".csv" )
  out.open_file
  taxonomy_levels = {}

  for i in 0...(names.size-1) do
    name_i = names[ i ]
    organism_i = @seqs[ name_i ].common_name

    for j in (i+1)...names.size do
      name_j = names[ j ]
      organism_j = @seqs[ name_j ].common_name

      out.write( seqs_file + "," )
      out.write( @seqs[ name_i ].sequence_name + "," + @seqs[ name_j ].sequence_name + "," )
      common_tax = taxons[ name_i ] & taxons[ name_j ]
      last_tax = common_tax[ common_tax.size-1 ]

      out.write( "#{last_tax},#{organism_i},#{organism_j}," )
      identities, mismatches, compared = compare_aligned( @aligned[ name_i ], @aligned[ name_j ] )
      percent = 0
      percent = (identities * 100) / compared if compared > 0
      fraction = 0
      fraction = (compared * 100) / @longest if @longest > 0
      # n = calculate_n( mismatches, v, compared, @q )
      n = calculate_n( mismatches, v, residues, @q )
      # q = calculate_q( mismatches, v, compared, n )
      q = calculate_q( mismatches, v, residues, n )
      out.write( "#{identities},#{mismatches},#{compared},#{percent}%,#{fraction}%,#{v},#{q},#{(n/2.0)},#{common_tax.join('|')}\n" )

      # Tally the results by common taxonomy for full length sequences.
      if ( compared + 100 >= max_length ) && ( percent >= MIN_PERCENT_IDENTITY )
        if taxonomy_levels[ last_tax ].nil?
          taxonomy_levels[ last_tax ] = TaxonomyLevel.new() 
          taxonomy_levels[ last_tax ].level = last_tax
          taxonomy_levels[ last_tax ].taxonomy = common_tax.join( "|" )
        end  # if

        taxonomy_levels[ last_tax ].add_comparison( identities, mismatches, compared, organism_i, organism_j )
      end  # if
    end  # for
  end  # for
  out.close_file

  taxonomy_levels.each do |name, tax_level|
    outliers = tax_level.flag_names
    puts "Outlier: #{seqs_file}.#{name}: #{outliers}" if outliers.size > 0
  end  # do

  # Estimate for Eutheria.
  level = "Eutheria"
  # level = "Mammalia"
  if (! taxonomy_levels[ level ].nil? )
    q_ave, q_max = calculate_q_range( taxonomy_levels[level], v, 200, 200 )
    # puts "Eutheria calculation 1: #{level}, N = 100, q = #{q_ave}"
    # puts "Eutheria calculation 2: #{level}, N = 100, q = #{q_max}"
    @q = q_max
    q_ave_eutheria = q_ave
  else
    return
  end  # if

  # Write out the Taxonomy summary information.
  out = OutputFile.new( seqs_file + "_tax.csv" )
  out.open_file
  taxonomy_levels.keys.each do |level|
    out.write( seqs_file + "," )
    out.write( taxonomy_levels[ level ].taxonomy + "," )
    out.write( level + ",#{v}," )
    out.write( "#{taxonomy_levels[ level ].percent}%," )
    residues = taxonomy_levels[level].residues
    mismatches = taxonomy_levels[level].mismatches
    n_ave, n_max = calculate_n_range( taxonomy_levels[level], v, q_ave_eutheria, @q )
    q_ave, q_max = calculate_q_range( taxonomy_levels[level], v, n_ave, n_max )
    names = taxonomy_levels[ level ].get_names.join( "|" )
    percent = 0
    percent = 100 - (mismatches * 100) /residues if residues > 0
    out.write( "#{mismatches}/#{residues},#{percent}%,#{n_ave/2},#{n_max/2},#{@variable.to_i},#{mismatches},#{residues},#{taxonomy_levels[level].names.keys.size},#{names},#{q_ave},#{q_max}\n" )
    puts "#{mismatches}/#{residues},#{percent}%,#{level} N_ave: #{(n_ave/2)} to #{(n_max/2)}; q #{q_ave} to #{q_max}"
    add_taxonomy_node( @taxonomy_tree, taxonomy_levels[ level ].taxonomy, level, n_ave, n_max, names, 0 )
  end  # do
  out.close_file

  eutherian_genes = 0
  eutherian_genes = taxonomy_levels[ "Eutheria" ].names.keys.size if ! taxonomy_levels[ "Eutheria" ].nil?
  write_v( seqs_file, v, q_ave_eutheria, first, eutherian_genes )

  puts "#{@taxonomy_tree}"
  puts "#{@taxonomy_leaves}"
end  # write_comparisons


################################################################################
def write_v( seqs_file, v, q_ave_eutheria, gene_name, eutherians )
  if ( v > 0.0 )
    out = OutputFile.new( seqs_file + "_V.csv" )
    out.open_file
    out.write( seqs_file + ",#{v},#{@q},#{q_ave_eutheria},#{count_cpg},#{@excluded_spacers},#{@chromosome[seqs_file]},#{@map[seqs_file]},#{eutherians}\n" )
    out.close_file
  end  # if
end  # write_v

################################################################################
def write_msa( seqs_file )
  names = @seqs.keys()
  first = names[0]
  return if @seqs[ first ].nil?
  first_seq = @seqs[ first ].sequence_data

  # Set up the alignment for the first sequence.
  @msa[first] = {}
  for i in 0...first_seq.length do
    @msa[first][i+1] = first_seq[i]
  end  # for 

  # Write out the alignment.
  out = OutputFile.new( seqs_file + ".fa" )
  out.open_file
  inserts = @tuples[ first ].inserts
  gap_widths = size_inserts( inserts )
  for name in names do
    @aligned[name] = ""

    if ! @msa[name].nil?
      align_seq = add_inserts( @msa[name], inserts, gap_widths, name )
      if ! align_seq.nil?
        out.write( ">" + @seqs[name].sequence_name + " " + @seqs[name].sequence_description + "\n" )
        @aligned[name] = align_seq.upcase
        out.write( "#{SeqTools::to_blocks( align_seq )}" )
      end  # if
    end  # if
  end  # for

  out.close_file
end  # write_msa


################################################################################
def print_msa
  names = @seqs.keys()
  first = names[0]
  first_seq = @seqs[ first ].sequence_data
  puts( first + "\t" + first_seq )
  for name in names[1..-1] do
    print( name + "\t" )
    for i in 1..@msa[name].keys.sort.last do
      if @msa[ name ][ i ].nil?
        print( "." )
      else
        print( @msa[ name ][ i ] )
      end  # if 
    end  # for
    print( "\n" )
  end  # for 

  # Data check:
  # for name in names[1..-1] do
  #   puts( name + "\t" + @msa[name].sort.join )
  # end  # for
  
  # Print out the insertions.
  puts
  inserts = @tuples[ first ].inserts
  inserts.keys.sort.each do |pos|
    inserts[pos].each do |name, ins|
      puts "#{pos}  #{name}  #{ins}"
    end  # do
  end  # do
end  # print_msa

################################################################################
def map_genes( mapfile )
  inf = InputFile.new( mapfile )
  inf.open_file()
  line = inf.next_line()   # header line
  while ( inf.is_end_of_file? == false )
    line = inf.next_line()
    if ( inf.is_end_of_file? == false )
      tokens = line.split
      gene_symbol = tokens[0].downcase
      @chromosome[ gene_symbol ] = tokens[1]
      # puts "#{gene_symbol}\t#{tokens[1]}"
      @map[ gene_symbol ] = tokens[2]
    end  # if
  end  # while
  inf.close_file()
end  # map_genes

################################################################################

end  # class Diverge


################################################################################
def test_msa( seqs_file )
  # puts "usage: ruby msa.rb msa2.fa\n" if ARGV.length < 1
  seqs_file = ARGV[0] if ARGV.length >= 1

  app = Diverge.new
  app.map_genes( "human.map" )
  seqs = app.load_fastas( seqs_file )		# protein FASTA sequences
  app.calculate_tuples
  app.compare_tuples 
  # app.print_msa
  app.write_msa( seqs_file )
  app.count_variable
  app.determine_cpg
  app.determine_spacers
  app.show_adaptation
  app.write_comparisons( seqs_file )
  app.count_residues( seqs_file )
  app.write_taxonomy_json( seqs_file )
  app.write_taxonomy_newick( seqs_file )
end  # method test_msa


################################################################################
  test_msa( "f9" )
# test_msa( "tp53" )
# test_msa( "serine2.fa" )
