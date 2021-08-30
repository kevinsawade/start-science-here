#!/usr/bin/env python3

import filecmp
test = filecmp.cmp('ingredients.txt', 'result.txt')

if not test:
    print("Something went wrong. Your ingredients.txt after the merge contains these lines:")
    with open('ingredients.txt', 'r') as f:
        print(f.read())
    print("\nBut the expected result should be:")
    with open('result.txt', 'r') as f:
        print(f.read())
    print("\nTry again!")
else:
    print("Success. You successfully carried out a 3-way merge and added more cilantro to the ingredients.txt")

