- N = 62500, d = 32, CPUs = 16, asynchronous

| n x t \ p | 2 | 4 | 6 | 8 | 12 | 16 |
| --- | --- |  --- | --- | --- | --- | --- |
| 1 x 16 | 70.7s | 35.1s | 24.48s | 20.63s | 14.64s | 12.35s |
| 2 x 8 | 72.6s | 33.44s | 27.51s | 20.85s | 14.79s | 12.24s |
| 4 x 4 | >73s | 35.6s | 27.38s | 21.98s | 16.28s | 12.35s |
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
