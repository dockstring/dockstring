"""Simple example of invoking dockstring."""

from dockstring import load_target

if __name__ == "__main__":
    target = load_target("DRD2")
    score, _ = target.dock("CC1=C(C(=O)N2CCCCC2=N1)CCN3CCC(CC3)C4=NOC5=C4C=CC(=C5)F")
    print(f"Docking completed without error, score={score:.3g}. END OF SCRIPT.")
