- N = 62500, d = 32, CPUs = 16, asynchronous

| n x t \ p | 2 | 4 | 6 | 8 | 12 | 16 |
| --- | --- |  --- | --- | --- | --- | --- |
| 1 x 16 | 70.7 | 35.1 | 24.48 | 20.63 | 14.64 | 12.35 |
| 2 x 8 | 72.6 | 33.44 | 27.51 | 20.85 | 14.79 | 12.24 |
| 4 x 4 | >73 | 35.6 | 27.38 | 21.98 | 16.28 | 12.35 |
| 8 x 2 | PD | PD | PD | PD | PD | PD |
| 16 x 1 | PD | PD | PD | PD | PD | PD |


- N = 62500, d = 32, CPUs = 16, synchronous

| n x t \ p | 2 | 4 | 6 | 8 | 12 | 16 |
| --- | --- |  --- | --- | --- | --- | --- |
| 1 x 16 | PD | PD | PD | PD | PD | PD |
| 2 x 8 | PD | PD | PD | PD | PD | PD |
| 4 x 4 | PD | PD | PD | PD | PD | PD |
| 8 x 2 | PD | PD | PD | PD | PD | PD |
| 16 x 1 | PD | PD | PD | PD | PD | PD |
