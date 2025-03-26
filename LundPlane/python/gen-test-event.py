#!/usr/bin/env python3
#
# Script to help produce simple small events for testing the
# azimuthal angle aspects of the LundPlane contrib
#
# Usage (from the LundPlane directory):
#
#   python3 python/gen-test-event.py | ./example_dpsi_collinear
#

import fastjet as fj
from math import *
import argparse

# produce an event with a collinear branching on each side,
# some angle phi between the two branches, and an additional
# particle in the z>0 hemisphere at phi!=0 (to test whether
# signs are set properly across the two hemispheres).
#
# (That particle generates the highest-kt primary on that side,
# but that primary does not pass the z-cut. Hence the two
# opposite-side primaries will have psibar != 0)
# 
# - smaller value of scale make the event more collinear
# - phi = 0 causes the soft particles to be on opposite sides
# - mirror_event rotates the event by pi (to test whether
#   seeding the declustering with the other jet makes a difference)

argparser = argparse.ArgumentParser()
argparser.add_argument("--scale", type=float, default=0.001)
argparser.add_argument("--phi", type=float, default=0.1)
argparser.add_argument("--mirror", action="store_true")
argparser.add_argument("--phi-ref", type=float, default=1.15)
args = argparser.parse_args()

scale = args.scale
phi   = args.phi
mirror_event = args.mirror


rap = lambda pt,z: -log(pt/z/2)

z1 = 0.12
pt1 = 0.007*scale
y1 = rap(pt1, z1)
phi1 = args.phi

# if args.phi_ref == 0, then we expect that
# psibar2 = phi2 - pi
z2 = 0.17
pt2 = 0.004*scale
y2 = -rap(pt2, z2)
phi2 = pi 

# highest kt emission, which will set psibar = 0
z3  = 0.05
pt3 = scale
y3  = rap(pt3, z3)
phi3 = args.phi_ref

particles_in = [
    fj.PseudoJet(0,0,+1,1),
    fj.PseudoJet(0,0,-1,1),
    fj.PtYPhiM(pt1, y1, phi1),
    # the psibar that will come out here should be 
    fj.PtYPhiM(pt2, y2, phi2),
    fj.PtYPhiM(pt3, y3, phi3),
]

# dumb way to get sum of all momenta -- clustering with e+e- Cambridge/Aachen
# with a large radius
psum = fj.JetDefinition(fj.ee_genkt_algorithm, 2*pi, 0.0)(particles_in)[0]
#print(psum)

# then boost to balance the event
for p in particles_in:
    p.unboost(psum)

# and mirror the event if requested
particles = []
for p in particles_in:
    if mirror_event:
        p = fj.PseudoJet(-p.px(), -p.py(), -p.pz(), p.E())
    particles.append(p)

for p in particles:
    print(p.px(), p.py(), p.pz(), p.E(), p.rap())
