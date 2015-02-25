# Bioviz-Sequence 

For this assignment I decided to align protein sequences. I chose an enzyme called *triphosphate nucleotidohydrolase*, which is involved in nucleotide metabolism. I found multiple versions of this enzyme using BLAST, each of which comes from a different organism. It could be interesting to see how the same enzyme is constructed differently in different organisms, especially when taking into account the phylogenetic connection between them.

#### Alignment
The alignment algorithm is an implementation of the algorithm form the tutorial included in the assignment (I made sure of this by running the sequences included in it). I then expanded it to align three sequences instead of two, which works well but led to a significant decrease in performance. Because the algorithm changed from O(m*n) to O(m*n*p), adding a third sequence of length 100 could slow it down by 100 times. While the three-sequence alignment is useful and interesting, it did have the downside of forcing me to use relatively short sequences (e.g. 150 characters).

#### Display
I decided to start with the most basic visualization of these sequences: displaying them as text stacked on top of each other. I then highlighted matches with blue and mismatches with red, to make it very apparent where each of these were (gaps were indicated by a '-' with no color, making them equally apparent). The use of color worked well because larger matched groups were made easily apparent by the contiguous areas of color.

To allow for arbitrarily long input, I made the text scale to fit the entire sequences. While this would eventually make the text unreadable, the relationship between the sequences would still be apparent from the colors. I added the ability to highlight a section of the sequences and see that in detail, which would allow the user to investigate the exact contents of an area of interest even in a very large dataset.

#### Free-Form
The tech section comprises the implementation of the following features:
* 3-sequence alignment.
* Sequence scaling for large input.
* Detail highlighting.
* Coloring based on matches.
 
The bio section comprises the usefulness of the following features in the comparison of  sequences such as the one included in the example (protein enzyme sequences):
* Ability to directly compare 3 sequences.
* Ability to see both the trends in alignment (matching/etc) and the specific characters that are aligned.
* Ability to investigate a section in detail.

#### Running It
This can be run by opening index.html with your browser. A live version is also available [here](http://jdb175.github.io/Bioviz-Sequence/). It might be a little slow to load, because of the algorithm performance issues outlined above.