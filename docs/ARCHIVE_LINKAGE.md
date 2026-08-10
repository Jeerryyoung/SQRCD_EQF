# Archive linkage and version policy

## Roles

1. **OSF registration** (`https://osf.io/n6gjy/`) records the prespecified plan and dated amendments. A frozen registration should not be rewritten to match later results.
2. **GitHub** (`https://github.com/Jeerryyoung/SQRCD_EQF`) is the living source-code repository.
3. **Zenodo** (`https://doi.org/10.5281/zenodo.21873862`) is the DOI-bearing archive for immutable research snapshots.

## Required reciprocal links

- GitHub README and `CITATION.cff` link to OSF and Zenodo.
- The Zenodo record should list the GitHub repository and OSF registration as related identifiers.
- The OSF project or registration metadata should list the GitHub repository and Zenodo DOI as related materials. If the registration is frozen, add the links to the associated OSF project or a dated amendment rather than altering the original registration answers.
- The manuscript Data and code availability statement should cite all three resources and the exact release tag.

## Submission freeze

1. Pass tests and remove sensitive or non-redistributable inputs.
2. Create a Git tag such as `v1.0.0-submission`.
3. Archive that exact tag on Zenodo and obtain its version-specific DOI.
4. Record the tag, commit SHA and DOI in the manuscript.
5. Keep the concept DOI for discovery and use the version DOI for the exact submitted snapshot.
