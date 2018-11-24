# @author Itrat Ahmed Akhter
# Implementation of functions and gradients needed by our examples. 
# They can return either interval or point evaluations depending on 
# whether the arguments are points or intervals
# This file also contains functions that returns linear constraints
# bounding the function given an interval bound. 
import math
import numpy as np
from scipy.spatial import ConvexHull
from intervalBasics import *


#tanh(a*x + b)
def tanhFun(x, a, b):
	tanhVal = np.tanh(a*x + b)
	if interval_p(x):
		lVal = np.nextafter(min(tanhVal[0], tanhVal[1]), float("-inf"))
		uVal = np.nextafter(max(tanhVal[0], tanhVal[1]), float("inf"))
		return np.array([lVal, uVal])

	return tanhVal

def tanhFunder(x, a, b):
	den = np.cosh(a*x + b)*np.cosh(a*x + b)
	grad = np.divide(a,den)
	separX = b/(-a*1.0)
	if interval_p(x):
		if (x[0] - separX)*(x[1] - separX) >= 0:
			grad = np.array([min(grad[0], grad[1]), max(grad[0], grad[1])])
		else:
			den0 = np.cosh(separX)*np.cosh(separX)
			grad0 = np.divide(a,den0)
			grad = np.array([min(grad[0], grad[1], grad0), max(grad[0], grad[1], grad0)])
		grad = interval_round(grad)

	return grad


#linear constraints in the form of a string for tanh
def tanhLinearConstraints(a, b, inputVar, outputVar, inputLow, inputHigh):
	separX = b/(-a*1.0)
	if a < 0:
		if inputLow <= separX and inputHigh <= separX:
			return triangleBounds(tanhFun, tanhFunder, inputVar, outputVar, inputLow, inputHigh, "neg", a, b)
		if inputLow >= separX and inputHigh >= separX:
			return triangleBounds(tanhFun, tanhFunder, inputVar, outputVar, inputLow, inputHigh, "pos", a, b)
	elif a >= 0:
		if inputLow <= separX and inputHigh <= separX:
			return triangleBounds(tanhFun, tanhFunder, inputVar, outputVar, inputLow, inputHigh, "pos", a, b)
		if inputLow >= separX and inputHigh >= separX:
			return triangleBounds(tanhFun, tanhFunder, inputVar, outputVar, inputLow, inputHigh, "neg", a, b)

	overallConstraint = "1 " + inputVar + " >= " + str(inputLow) + "\n"
	overallConstraint += "1 " + inputVar + " <= " + str(inputHigh) + "\n"
	overallConstraint += "1 " + outputVar + " <= 1.0\n"
	overallConstraint += "1 " + outputVar + " >= -1.0\n"
	allTrianglePoints = []
	allTrianglePoints += trianglePoints(tanhFun, tanhFunder, inputLow, 0.0, a, b)
	allTrianglePoints += trianglePoints(tanhFun, tanhFunder, 0.0, inputHigh, a, b)
	allTrianglePoints = np.array(allTrianglePoints)
	try:
		cHullConstraints = convexHullConstraints2D(allTrianglePoints, inputVar, outputVar)
		overallConstraint += cHullConstraints
	except:
		pass
	return overallConstraint


# This function calculates tangents at inputLow and inputHigh and
# and finds the intersection between the tangents. 
# It returns the three points of a triangle:
# (inputLow, function(inputLow)), (inputHigh, function(inputHigh)), (intersectionX, function(intersectionX))
def trianglePoints(function, functionDer, inputLow, inputHigh, a, b=None):
	if inputLow > inputHigh:
		return []
	if b is None:
		funLow = function(inputLow, a)
		dLow = functionDer(inputLow, a)
		funHigh = function(inputHigh, a)
		dHigh = functionDer(inputHigh, a)
	else:
		funLow = function(inputLow, a, b)
		dLow = functionDer(inputLow, a, b)
		funHigh = function(inputHigh, a, b)
		dHigh = functionDer(inputHigh, a, b)
	
	cLow = funLow - dLow*inputLow
	cHigh = funHigh - dHigh*inputHigh

	diff = inputHigh - inputLow
	if(diff == 0):
		diff = 1e-10
	dThird = (funHigh - funLow)/diff
	cThird = funLow - dThird*inputLow

	leftIntersectX, leftIntersectY = None, None
	if abs(dHigh - dLow) < 1e-8:
		leftIntersectX = inputLow
		leftIntersectY = funLow
	else:
		leftIntersectX = (cHigh - cLow)/(dLow - dHigh)
		leftIntersectY = dLow*leftIntersectX + cLow

	#print ("leftIntersectX", leftIntersectX, "leftIntersectY", leftIntersectY)
	tPts = [[inputLow, funLow],[inputHigh, funHigh],[leftIntersectX, leftIntersectY]]
	return tPts


# This function constructs linear constraints from the given interval bounds
# If there are no inflection points of the function in the interval bounds
# then it returns triangle constraints formed from the tangents at the interval bounds
# and a secant line between the interval bounds depending on the convexity of the
# function. Otherwise it just returns constraints indicating the interval bounds
def triangleBounds(function, functionDer, inputVar, outputVar, inputLow, inputHigh, secDer, a, b=None):
	if b is None:
		funLow = function(np.array([inputLow]), a)[0]
		dLow = functionDer(np.array([inputLow]), a)[0]
		funHigh = function(np.array([inputHigh]), a)[0]
		dHigh = functionDer(np.array([inputHigh]), a)[0]
	else:
		funLow = function(np.array([inputLow]), a, b)[0]
		dLow = functionDer(np.array([inputLow]), a, b)[0]
		funHigh = function(np.array([inputHigh]), a, b)[0]
		dHigh = functionDer(np.array([inputHigh]), a, b)[0]
	
	cLow = funLow - dLow*inputLow
	cHigh = funHigh - dHigh*inputHigh

	diff = inputHigh - inputLow
	if(diff == 0):
		diff = 1e-10
	dThird = (funHigh - funLow)/diff
	cThird = funLow - dThird*inputLow

	overallConstraint = ""
	overallConstraint += "1 " + inputVar + " >= " + str(inputLow) + "\n"
	overallConstraint += "1 " + inputVar + " <= " + str(inputHigh) + "\n"
	if secDer == None:
		return overallConstraint

	if secDer == "pos":
		return overallConstraint + "1 "+ outputVar + " + " +str(-dThird) + " " + inputVar + " <= "+str(cThird)+"\n" +\
				"1 "+outputVar + " + " +str(-dLow) + " " + inputVar + " >= "+str(cLow)+"\n" +\
				"1 "+outputVar + " + " +str(-dHigh) + " " + inputVar + " >= "+str(cHigh) + "\n"

	
	if secDer == "neg":
		return overallConstraint + "1 "+ outputVar + " + " +str(-dThird) + " " + inputVar + " >= "+str(cThird)+"\n" +\
				"1 "+outputVar + " + " +str(-dLow) + " " + inputVar + " <= "+str(cLow)+"\n" +\
				"1 "+outputVar + " + " +str(-dHigh) + " " + inputVar + " <= "+str(cHigh) + "\n"


# This function finds the convex hull of a list of 2d points and creates
# constraints around the convex hull and returns it in the form of strings
def convexHullConstraints2D(points, inputVar, outputVar):
	#print ("points")
	#print (points)
	hull = ConvexHull(points)
	convexHullMiddle = np.zeros((2))
	numPoints = 0
	for simplex in hull.simplices:
		#print ("simplex", simplex)
		for ind in simplex:
			convexHullMiddle += [points[ind,0],points[ind,1]]
			numPoints += 1
	convexHullMiddle = convexHullMiddle/(numPoints*1.0)
	#print ("convexHullMiddle", convexHullMiddle)
	overallConstraint = ""
	for si in range(len(hull.simplices)):
		simplex = hull.simplices[si]
		#print ("simplex", simplex)

		pt1x = points[simplex[0],0]
		pt1y = points[simplex[0],1]

		pt2x = points[simplex[1],0]
		pt2y = points[simplex[1],1]

		#print ("pt1x ", pt1x, "pt1y", pt1y)
		#print ("pt2x ", pt2x, "pt2y", pt2y)

		grad = (pt2y - pt1y)/(pt2x - pt1x)
		c = pt1y - grad*pt1x
		#print ("grad", grad, "c", c)

		yMiddle = grad*convexHullMiddle[0] + c
		#print ("yMiddle", yMiddle)

		sign = " <= "
		if convexHullMiddle[1] > yMiddle:
			sign = " >= "

		#print ("sign", sign)

		overallConstraint += "1 " + outputVar + " + " + str(-grad) + " " + inputVar + \
			sign + str(c) + "\n"
		
	#print ("overallConstraint")
	#print (overallConstraint)
	return overallConstraint
