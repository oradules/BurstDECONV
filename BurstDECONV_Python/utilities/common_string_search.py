#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 10:31:08 2022

@author: mdouaihy
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 10:06:01 2022

@author: mdouaihy
"""


# Python code for the above approach
from functools import lru_cache
from operator import itemgetter

def longest_common_substring(x: str, y: str) -> (int, int, int):
	
	# function to find the longest common substring

	# Memorizing with maximum size of the memory as 1
	@lru_cache(maxsize=1)
	
	# function to find the longest common prefix
	def longest_common_prefix(i: int, j: int) -> int:
	
		if 0 <= i < len(x) and 0 <= j < len(y) and x[i] == y[j]:
			return 1 + longest_common_prefix(i + 1, j + 1)
		else:
			return 0

	# diagonally computing the subproblems
	# to decrease memory dependency
	def digonal_computation():
		
		# upper right triangle of the 2D array
		for k in range(len(x)):	
			yield from ((longest_common_prefix(i, j), i, j)
						for i, j in zip(range(k, -1, -1),
									range(len(y) - 1, -1, -1)))
		
		# lower left triangle of the 2D array
		for k in range(len(y)):	
			yield from ((longest_common_prefix(i, j), i, j)
						for i, j in zip(range(k, -1, -1),
									range(len(x) - 1, -1, -1)))

	# returning the maximum of all the subproblems
	return max(digonal_computation(), key=itemgetter(0), default=(0, 0, 0))


def common_substring_in_list(llist):
    base_word = llist[0].replace('CalibratedTraces','')
    base_word = base_word.replace('.npz', '')
    base_word = base_word.replace('result_', '')
    base_word = base_word.lower()

    for i in range(1,len(llist)):
        new_world = llist[i].lower()
        length, i, j = longest_common_substring(base_word, new_world)
        base_word = base_word[i: i + length]
    return llist[-1][j: j+length]

def common_first_two_substrings_in_list(llist):
    first_base_common = common_substring_in_list(llist)
    first_common = first_base_common.lower()
    
    new_list = []
    for i in range(len(llist)):
        new_word = llist[i].replace('CalibratedTraces', '')
        new_word = new_word.lower().replace(first_common.lower(), '')
        
        new_list.append(new_word)
    
    second_common = common_substring_in_list(new_list)
    return(second_common + '_' + first_base_common)