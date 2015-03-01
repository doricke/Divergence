
require 'fasta_sequence.rb'
require 'fasta_iterator.rb'
require 'output_file.rb'

################################################################################
# Copyright (C) 2015  Darrell O. Ricke, PhD
# Author::	Darrell O. Ricke, Ph.D.	(mailto: Darrell.Ricke@ll.mit.edu)
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

###############################################################################
class SplitByGene

###############################################################################

attr_accessor :genes  	    # genes 
attr_accessor :genes_human  # human genes


###############################################################################
def initialize
  @genes = {}
  @genes_human = {}
end  # method initialize


###############################################################################
def add_gene( fasta )
  taxonomy = fasta.annotation[ "taxonomy" ]
  gene_name = fasta.annotation[ "gene" ]

  # Select vertebrate organisms.
  # if ( taxonomy.include?( "Vertebrata" ) && ( ! gene_name.nil?) &&  ( gene_name.length > 0 ) )
  if ( ( ! gene_name.nil?) &&  ( gene_name.length > 0 ) )
    gene_name.delete!( "/,()'<>:;`@#\$%^*+={}[]?/\|,\." )
    gene_name.downcase!
    @genes[ gene_name ] = [] if @genes[ gene_name ].nil?

    organism = fasta.annotation[ "organism" ]
    if ( ! organism.nil? ) && ( organism.include?( "Homo sapiens" ) )
      @genes_human[ gene_name ] = fasta 
    else
      @genes[ gene_name ] << fasta 
    end  # if

    # Matching gene names for this organism.
    # resolve( gene1, fasta, gene1.rna_seq, gene2.rna_seq )
  end  # if
end  # method add_gene


###############################################################################
def read_sequences( filename )
  in_fasta = FastaIterator.new( filename )
  in_fasta.open_file
  while ( in_fasta.is_end_of_file? == false )
    fasta = in_fasta.next_sequence
    fasta.parse_annotation
    add_gene( fasta ) if ( ! fasta.nil? ) && ( fasta.sequence_data.length > 0 )
  end  # while

  in_fasta.close_file
end  # method read_sequences


###############################################################################
def resolve( gene1, gene2, words1, words2 )
  len = words1.size
  len = words2.size if words2.size < words1.size

  @organism_genes[ gene2.organism ][ gene2.name ] = gene2 if words2.size > words1.size 
end  # method resolve


###############################################################################
def write_genes
  @genes.keys.each do |gene_name|
    # if ( ! @genes_human[ gene_name ].nil? )
      filename = gene_name.gsub( " ", "_" )
      out = OutputFile.new( filename )
      out.open_file
      # out.write( @genes_human[ gene_name ].to_string )

      @genes[ gene_name ].each do |fasta|
        out.write( fasta.to_string )
      end  # do
      out.close_file
    # end  # if
  end  # do
end  # method write_genes


###############################################################################

end  # class SplitByGene


###############################################################################
def main_select_genes
  filename = "sprot"
  app = SplitByGene.new
  app.read_sequences( filename )
  app.write_genes
end  # method main_select_genes


###############################################################################

main_select_genes
