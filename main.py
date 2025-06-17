#!/usr/bin/env python3
import argparse
import numpy as np
from tssltnwyckoff import builda2dcrystal, generallatticevectors
from WyckoffPlaneGroupID import twodstructureplotter

def main():
    parser = argparse.ArgumentParser(
        description='Plot 2D crystal structures from plane groups and Wyckoff positions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py p1 a --size 3
  python main.py p2 a b --lattice --bonds
  python main.py p4m a b c --species --bonds --size 4

Available plane groups:
  p1, p2, pm, pg, cm, pmm, pmg, pgg, cmm, p4, p4m, p4g, p3, p3m, p3m1, p6, p6m
        """
    )

    # Required arguments
    parser.add_argument('planegroup', type=str, help='Plane group symbol (e.g., p1, p2, pm)')
    parser.add_argument('wyckoff_positions', type=str, nargs='+', help='Occupied Wyckoff positions (e.g., a b c)')

    # Optional visualization arguments
    parser.add_argument('--size', type=int, default=2, help='Size of the supercell (default: 2)')
    parser.add_argument('--lattice', action='store_true', help='Show lattice points')
    parser.add_argument('--species', action='store_true', help='Show species points')
    parser.add_argument('--bonds', action='store_true', help='Show bonds')
    parser.add_argument('--planegroup_options', action='store_true', help='List available Wyckoff positions for the given plane group and exit without plotting')

    args = parser.parse_args()

    # If planegroup_options is set, print available Wyckoff positions and exit
    if args.planegroup_options:
        print(f"Available Wyckoff positions for plane group {args.planegroup}:")
        
        wyckoff_positions = {
            "p1": {
                "a": "[(0.5, 0.5)]"
            },
            "p2": {
                "a": "[(0, 0)]",
                "b": "[(0, 0.5)]",
                "c": "[(0.5, 0)]",
                "d": "[(0.5, 0.5)]",
                "e": "[(0.25, 0.25), (-0.25, -0.25)]"
            },
            "pm": {
                "a": "[(0, 0.5)]",
                "b": "[(0.5, 0.5)]",
                "c": "[(0.25, 0.25), (-0.25, 0.25)]"
            },
            "pg": {
                "a": "[(0.25, 0.25), (-0.25, 0.75)]"
            },
            "cm": {
                "a": "[(0, 0.5)]",
                "b": "[(0.25, 0.25), (-0.25, 0.25)]"
            },
            "pmm": {
                "a": "[(0, 0)]",
                "b": "[(0, 0.5)]",
                "c": "[(0.5, 0)]",
                "d": "[(0.5, 0.5)]",
                "e": "[(0.75, 0), (-0.75, 0)]",
                "f": "[(0.75, 0.5), (-0.75, 0.5)]",
                "g": "[(0, 0.75), (0, -0.75)]",
                "h": "[(0.5, 0.75), (0.5, -0.75)]",
                "i": "[(0.25, 0.25), (-0.25, -0.25), (-0.25, 0.25), (0.25, -0.25)]"
            },
            "pmg": {
                "a": "[(0, 0), (0.5, 0)]",
                "b": "[(0, 0.5), (0.5, 0.5)]",
                "c": "[(0.25, 0.25), (0.75, -0.25)]",
                "d": "[(0.125, 0.125), (0.375, 0.125), (-0.125, -0.125), (0.625, -0.125)]"
            },
            "pgg": {
                "a": "[(0, 0), (0.5, 0.5)]",
                "b": "[(0, 0.5), (0.5, 0)]",
                "c": "[(0.125, 0.125), (-0.125, -0.125), (0.375, 0.625), (0.625, 0.375)]"
            },
            "cmm": {
                "a": "[(0, 0), (0.5, 0.5)]",
                "b": "[(0, 0.5), (0.5, 1)]",
                "c": "[(0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (1.25, 0.75)]",
                "d": "[(0.25, 0), (-0.25, 0), (-0.75, 0), (0.75, 0)]",
                "e": "[(0, 0.25), (0, -0.25), (0, 0.75), (0, -0.75)]",
                "f": "[(0.125, 0.125), (-0.125, -0.125), (-0.125, 0.125), (0.125, -0.125), (0.625, 0.625), (-0.625, -0.625), (0.625, -0.625), (-0.625, 0.625)]"
            },
            "p4": {
                "a": "[(0, 0)]",
                "b": "[(0.5, 0.5)]",
                "c": "[(0.5, 0), (0, 0.5)]",
                "d": "[(0.25, 0.25), (-0.25, 0.25), (-0.25, -0.25), (0.25, -0.25)]"
            },
            "p4m": {
                "a": "[(0, 0)]",
                "b": "[(0.5, 0.5)]",
                "c": "[(0.5, 0), (0, 0.5)]",
                "d": "[(0.29, 0), (-0.29, 0), (0, -0.29), (0, 0.29)]",
                "e": "[(0.25, 0.5), (-0.25, 0.5), (0.5, 0.25), (0.5, -0.25)]",
                "f": "[(0.25, 0.25), (-0.25, 0.25), (-0.25, -0.25), (0.25, -0.25)]",
                "g": "[(0.25, 0.75), (-0.25, 0.75), (-0.25, -0.75), (0.25, -0.75), (0.75, 0.25), (-0.75, 0.25), (0.75, -0.25), (-0.75, -0.25)]"
            },
            "p4g": {
                "a": "[(0, 0), (0.5, 0.5)]",
                "b": "[(0, 0.5), (0.5, 0)]",
                "c": "[(0.25, 0.75), (-0.25, 0.25), (0.25, 0.25), (0.75, -0.25)]",
                "d": "[(0.125, 0.75), (0.75, -0.125), (-0.125, -0.75), (-0.75, 0.125), (0.375, 1.25), (0.625, -0.25), (1.25, 0.625), (-0.25, -0.375)]"
            },
            "p3": {
                "a": "[(0, 0)]",
                "b": "[(0.3333, 0.6667)]",
                "c": "[(0.6667, 0.3333)]",
                "d": "[(0.125, 0.75), (-0.75, -0.625), (0.625, -0.125)]"
            },
            "p3m": {
                "a": "[(0, 0)]",
                "b": "[(0.3333, 0.6667)]",
                "c": "[(0.6667, 0.3333)]",
                "d": "[(0.25, -0.25), (0.25, 0.5), (-0.5, -0.25)]",
                "e": "[(0.125, 0.75), (-0.75, -0.625), (0.625, -0.125), (-0.75, -0.125), (0.625, 0.75), (0.125, -0.625)]"
            },
            "p3m1": {
                "a": "[(0, 0)]",
                "b": "[(0.3333, 0.6667), (0.6667, 0.3333)]",
                "c": "[(0.25, 0), (0, 0.25), (-0.25, -0.25)]",
                "d": "[(0.125, 0.75), (-0.75, -0.625), (0.625, -0.125), (0.75, 0.125), (-0.625, -0.75), (-0.125, 0.625)]"
            },
            "p6": {
                "a": "[(0, 0)]",
                "b": "[(0.3333, 0.6667), (0.6667, 0.3333)]",
                "c": "[(0.5, 0), (0, 0.5), (0.5, 0.5)]",
                "d": "[(0.125, 0.75), (-0.75, -0.625), (0.625, -0.125), (-0.125, -0.75), (0.75, 0.625), (-0.625, 0.125)]"
            },
            "p6m": {
                "a": "[(0, 0)]",
                "b": "[(0.3333, 0.6667), (0.6667, 0.3333)]",
                "c": "[(0.5, 0), (0, 0.5), (0.5, 0.5)]",
                "d": "[(0.375, 0), (0, 0.375), (-0.375, -0.375), (-0.375, 0), (0, -0.375), (0.375, 0.375)]",
                "e": "[(0.25, -0.25), (0.25, 0.5), (-0.5, -0.25), (-0.25, 0.25), (-0.25, -0.5), (0.5, 0.25)]",
                "f": "[(0.125, 0.75), (-0.75, -0.625), (0.625, -0.125), (-0.125, -0.75), (0.75, 0.625), (-0.625, 0.125), (-0.75, -0.125), (0.625, 0.75), (0.125, -0.625), (0.75, 0.125), (-0.625, -0.75), (-0.125, 0.625)]"
            }
        }
        
        if args.planegroup in wyckoff_positions:
            for letter, coords in wyckoff_positions[args.planegroup].items():
                print(f"  {letter}: {coords}")
        else:
            print(f"Error: Plane group '{args.planegroup}' not found.")
            print("Available plane groups:")
            print("  " + ", ".join(wyckoff_positions.keys()))
        
        exit(0)

    # Build crystal structure
    fraccoords = builda2dcrystal(args.planegroup, args.wyckoff_positions)
    print(f"Generated fractional coordinates:\n{fraccoords}")
    
    # Get lattice vectors for the plane group
    lattice_vectors = generallatticevectors(args.planegroup)
    
    # Convert fractional coordinates to numpy array
    cartpts = np.zeros_like(fraccoords)#turns form fractional into cartesian coords
    for i in range(fraccoords.shape[0]):
        cartpts[i,:] = lattice_vectors[0,:2]*fraccoords[i,0] + lattice_vectors[1,:2]*fraccoords[i,1]
    
    # Plot the structure with specified options
    twodstructureplotter(
        lattice_vectors, 
        cartpts, 
        args.size,
        latticeshown=args.lattice,
        speciesshown=args.species,
        bondshown=args.bonds
    )

if __name__ == '__main__':
    main()
