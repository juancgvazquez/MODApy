import pandas as pa


def rmerge(left, right, **kwargs):
	"""Perform a merge using pandas with optional removal of overlapping
	column names not associated with the join.

	Though I suspect this does not adhere to the spirit of pandas merge
	command, I find it useful because re-executing IPython notebook cells
	containing a merge command does not result in the replacement of existing
	columns if the name of the resulting DataFrame is the same as one of the
	two merged DataFrames, i.e. data = pa.merge(data,new_dataframe). I prefer
	this command over pandas df.combine_first() method because it has more
	flexible join options.

	The column removal is controlled by the 'replace' flag which is
	'left' (default) or 'right' to remove overlapping columns in either the
	left or right DataFrame. If 'replace' is set to None, the default
	pandas behavior will be used. All other parameters are the same
	as pandas merge command.

	Examples
	--------
	>>> left       >>> right
	   a  b   c       a  c   d
	0  1  4   9    0  1  7  13
	1  2  5  10    1  2  8  14
	2  3  6  11    2  3  9  15
	3  4  7  12

	>>> rmerge(left,right,on='a')
	   a  b  c   d
	0  1  4  7  13
	1  2  5  8  14
	2  3  6  9  15

	>>> rmerge(left,right,on='a',how='left')
	   a  b   c   d
	0  1  4   7  13
	1  2  5   8  14
	2  3  6   9  15
	3  4  7 NaN NaN

	>>> rmerge(left,right,on='a',how='left',replace='right')
	   a  b   c   d
	0  1  4   9  13
	1  2  5  10  14
	2  3  6  11  15
	3  4  7  12 NaN

	>>> rmerge(left,right,on='a',how='left',replace=None)
	   a  b  c_x  c_y   d
	0  1  4    9    7  13
	1  2  5   10    8  14
	2  3  6   11    9  15
	3  4  7   12  NaN NaN
	"""

	# Function to flatten lists from http://rosettacode.org/wiki/Flatten_a_list#Python
	def flatten(lst):
		return sum(([x] if not isinstance(x, list) else flatten(x) for x in lst), [])

	# Set default for removing overlapping columns in "left" to be true
	myargs = {'replace': 'left'}
	myargs.update(kwargs)

	# Remove the replace key from the argument dict to be sent to
	# pandas merge command
	kwargs = {k: v for k, v in myargs.items() if k is not 'replace'}

	if myargs['replace'] is not None:
		# Generate a list of overlapping column names not associated with the join
		skipcols = set(flatten([v for k, v in myargs.items() if k in ['on', 'left_on', 'right_on']]))
		leftcols = set(left.columns)
		rightcols = set(right.columns)
		dropcols = list((leftcols & rightcols).difference(skipcols))

		# Remove the overlapping column names from the appropriate DataFrame
		if myargs['replace'].lower() == 'left':
			left = left.copy().drop(dropcols, axis=1)
		elif myargs['replace'].lower() == 'right':
			right = right.copy().drop(dropcols, axis=1)

	df = pa.merge(left, right, **kwargs)

	return df