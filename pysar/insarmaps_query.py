#! /usr/bin/env python

import urllib2
import getopt
import os
import argparse

class BasicHTTP:
	@staticmethod
	def get(url):
		res = urllib2.urlopen(url)
		return res.read()

def buildURL(args):
	url = "http://ec2-52-41-231-16.us-west-2.compute.amazonaws.com/WebServices?"

	if args.dataset:
		url += "dataset=" + args.dataset + "&"
	if args.longitude:
		url += "longitude=" + args.longitude + "&"
	if args.latitude:
		url += "latitude=" + args.latitude + "&"

	return url[:-1]

def main():
	parser = argparse.ArgumentParser(description='Query insarmaps database.')
	parser.add_argument("-s", "--satellite", help="satellite to search for")
	parser.add_argument("-r", "--relativeOrbit", help="relative orbit to search for")
	parser.add_argument("-f", "--firstFrame", help="first frame to search for")
	parser.add_argument("-m", "--mode", help="mode to search for")
	parser.add_argument("-d", "--flightDir", help="flight direction to search for")
	parser.add_argument("-D", "--dataset", help="dataset to search in")
	parser.add_argument("-l", "--latitude", help="latitude of point to search for")
	parser.add_argument("-L", "--longitude", help="longitude of point to search for")

	parseArgs = parser.parse_args()

	url = buildURL(parseArgs)

	print BasicHTTP.get(url)

if __name__ == '__main__':
	main()
