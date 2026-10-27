"""Generate revised scientific figures and their supporting numerical data.

The original v1 figure generator is retained in baseline/recovered_v1_source.zip.
This entry point cannot produce old hard-coded artifact labels for new data.
"""
from .phase1 import run_study


def main():
    run_study()


if __name__=='__main__':main()
