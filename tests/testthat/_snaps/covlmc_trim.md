# trimmed model can be drawn

    Code
      draw(t_m_cov_model, model = "coef", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ (1.28,7.54]/(0,1.28] | -3 | -0.43 | -0.64 | 1.5 | -2.3 | 2.6 ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ (1.28,7.54]/(0,1.28] | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ (1.28,7.54]/(0,1.28] | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ (1.28,7.54]/(0,1.28] | -2.7 | 1.3 ])
      |   '-- (1.28,7.54] (0.15 [ (1.28,7.54]/(0,1.28] | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ (1.28,7.54]/(0,1.28] | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ (1.28,7.54]/(0,1.28] | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ (1.28,7.54]/(0,1.28] | 0.99 | 1.3 ])
                  '-- (1.28,7.54] (0.72 [ (1.28,7.54]/(0,1.28] | 2.8 ])

---

    Code
      draw(t_m_cov, model = "coef", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ (1.28,7.54]/(0,1.28] | -3 | -0.43 | -0.64 | 1.5 | -2.3 | 2.6 ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ (1.28,7.54]/(0,1.28] | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ (1.28,7.54]/(0,1.28] | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ (1.28,7.54]/(0,1.28] | -2.7 | 1.3 ])
      |   '-- (1.28,7.54] (0.15 [ (1.28,7.54]/(0,1.28] | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ (1.28,7.54]/(0,1.28] | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ (1.28,7.54]/(0,1.28] | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ (1.28,7.54]/(0,1.28] | 0.99 | 1.3 ])
                  '-- (1.28,7.54] (0.72 [ (1.28,7.54]/(0,1.28] | 2.8 ])

---

    Code
      draw(t_m_cov_model, model = "full", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ ((0,1.28])  | (I) | day_night_1TRUE | day_night_2TRUE | day_night_3TRUE | day_night_4TRUE | day_night_5TRUE
      |   |   |   |   |                      (1.28,7.54] | -3  | -0.43           | -0.64           | 1.5             | -2.3            | 2.6             ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ ((0,1.28])  | (I) 
      |   |   |   |                           (1.28,7.54] | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ ((0,1.28])  | (I) 
      |   |   |                          (1.28,7.54] | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ ((0,1.28])  | (I)  | day_night_1TRUE
      |   |                            (1.28,7.54] | -2.7 | 1.3             ])
      |   '-- (1.28,7.54] (0.15 [ ((0,1.28])  | (I) 
      |                           (1.28,7.54] | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ ((0,1.28])  | (I)
          |                    (1.28,7.54] | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ ((0,1.28])  | (I)
              |                 (1.28,7.54] | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ ((0,1.28])  | (I)  | day_night_1TRUE
                  |                     (1.28,7.54] | 0.99 | 1.3             ])
                  '-- (1.28,7.54] (0.72 [ ((0,1.28])  | (I)
                                          (1.28,7.54] | 2.8 ])

---

    Code
      draw(t_m_cov, model = "full", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ ((0,1.28])  | (I) | day_night_1TRUE | day_night_2TRUE | day_night_3TRUE | day_night_4TRUE | day_night_5TRUE
      |   |   |   |   |                      (1.28,7.54] | -3  | -0.43           | -0.64           | 1.5             | -2.3            | 2.6             ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ ((0,1.28])  | (I) 
      |   |   |   |                           (1.28,7.54] | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ ((0,1.28])  | (I) 
      |   |   |                          (1.28,7.54] | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ ((0,1.28])  | (I)  | day_night_1TRUE
      |   |                            (1.28,7.54] | -2.7 | 1.3             ])
      |   '-- (1.28,7.54] (0.15 [ ((0,1.28])  | (I) 
      |                           (1.28,7.54] | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ ((0,1.28])  | (I)
          |                    (1.28,7.54] | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ ((0,1.28])  | (I)
              |                 (1.28,7.54] | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ ((0,1.28])  | (I)  | day_night_1TRUE
                  |                     (1.28,7.54] | 0.99 | 1.3             ])
                  '-- (1.28,7.54] (0.72 [ ((0,1.28])  | (I)
                                          (1.28,7.54] | 2.8 ])

---

    Code
      draw(t_m_cov_model, model = "coef", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ (0.48,1.63]/(0,0.48] | -2.7 | 1.1 
      |                       (1.63,7.54]/(0,0.48] | -4.5 | 0.84 ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ (0.48,1.63]/(0,0.48] | 0.12  | -16 | 17
      |   |                     (1.63,7.54]/(0,0.48] | -0.37 | -16 | 15 ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ (0.48,1.63]/(0,0.48] | 2     | 0.58
      |                                          (1.63,7.54]/(0,0.48] | -0.31 | 1    ])
      '-- (1.63,7.54] (2.4e-06 [ (0.48,1.63]/(0,0.48] | 18 | -15
                                 (1.63,7.54]/(0,0.48] | 20 | -16 ])

---

    Code
      draw(t_m_cov, model = "coef", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ (0.48,1.63]/(0,0.48] | -2.7 | 1.1 
      |                       (1.63,7.54]/(0,0.48] | -4.5 | 0.84 ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ (0.48,1.63]/(0,0.48] | 0.12  | -16 | 17
      |   |                     (1.63,7.54]/(0,0.48] | -0.37 | -16 | 15 ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ (0.48,1.63]/(0,0.48] | 2     | 0.58
      |                                          (1.63,7.54]/(0,0.48] | -0.31 | 1    ])
      '-- (1.63,7.54] (2.4e-06 [ (0.48,1.63]/(0,0.48] | 18 | -15
                                 (1.63,7.54]/(0,0.48] | 20 | -16 ])

---

    Code
      draw(t_m_cov_model, model = "full", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ ((0,0.48])  | (I)  | day_night_1TRUE
      |                       (0.48,1.63] | -2.7 | 1.1            
      |                       (1.63,7.54] | -4.5 | 0.84            ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ ((0,0.48])  | (I)   | day_night_1TRUE | day_night_2TRUE
      |   |                     (0.48,1.63] | 0.12  | -16             | 17             
      |   |                     (1.63,7.54] | -0.37 | -16             | 15              ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ ((0,0.48])  | (I)   | day_night_1TRUE
      |                                          (0.48,1.63] | 2     | 0.58           
      |                                          (1.63,7.54] | -0.31 | 1               ])
      '-- (1.63,7.54] (2.4e-06 [ ((0,0.48])  | (I) | day_night_1TRUE
                                 (0.48,1.63] | 18  | -15            
                                 (1.63,7.54] | 20  | -16             ])

---

    Code
      draw(t_m_cov, model = "full", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ ((0,0.48])  | (I)  | day_night_1TRUE
      |                       (0.48,1.63] | -2.7 | 1.1            
      |                       (1.63,7.54] | -4.5 | 0.84            ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ ((0,0.48])  | (I)   | day_night_1TRUE | day_night_2TRUE
      |   |                     (0.48,1.63] | 0.12  | -16             | 17             
      |   |                     (1.63,7.54] | -0.37 | -16             | 15              ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ ((0,0.48])  | (I)   | day_night_1TRUE
      |                                          (0.48,1.63] | 2     | 0.58           
      |                                          (1.63,7.54] | -0.31 | 1               ])
      '-- (1.63,7.54] (2.4e-06 [ ((0,0.48])  | (I) | day_night_1TRUE
                                 (0.48,1.63] | 18  | -15            
                                 (1.63,7.54] | 20  | -16             ])

---

    Code
      draw(t_m_cov_model, model = "coef", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ 1/0 | -3 | -0.43 | -0.64 | 1.5 | -2.3 | 2.6 ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ 1/0 | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ 1/0 | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ 1/0 | -2.7 | 1.3 ])
      |   '-- (1.28,7.54] (0.15 [ 1/0 | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ 1/0 | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ 1/0 | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ 1/0 | 0.99 | 1.3 ])
                  '-- (1.28,7.54] (0.72 [ 1/0 | 2.8 ])

---

    Code
      draw(t_m_cov, model = "coef", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ 1/0 | -3 | -0.43 | -0.64 | 1.5 | -2.3 | 2.6 ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ 1/0 | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ 1/0 | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ 1/0 | -2.7 | 1.3 ])
      |   '-- (1.28,7.54] (0.15 [ 1/0 | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ 1/0 | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ 1/0 | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ 1/0 | 0.99 | 1.3 ])
                  '-- (1.28,7.54] (0.72 [ 1/0 | 2.8 ])

---

    Code
      draw(t_m_cov_model, model = "full", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ (0) | (I) | day_night_1TRUE | day_night_2TRUE | day_night_3TRUE | day_night_4TRUE | day_night_5TRUE
      |   |   |   |   |                      1   | -3  | -0.43           | -0.64           | 1.5             | -2.3            | 2.6             ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ (0) | (I) 
      |   |   |   |                           1   | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ (0) | (I) 
      |   |   |                          1   | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ (0) | (I)  | day_night_1TRUE
      |   |                            1   | -2.7 | 1.3             ])
      |   '-- (1.28,7.54] (0.15 [ (0) | (I) 
      |                           1   | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ (0) | (I)
          |                    1   | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ (0) | (I)
              |                 1   | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ (0) | (I)  | day_night_1TRUE
                  |                     1   | 0.99 | 1.3             ])
                  '-- (1.28,7.54] (0.72 [ (0) | (I)
                                          1   | 2.8 ])

---

    Code
      draw(t_m_cov, model = "full", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,1.28]
      |   +-- (0,1.28]
      |   |   +-- (0,1.28]
      |   |   |   +-- (0,1.28]
      |   |   |   |   +-- (0,1.28] (0.0025 [ (0) | (I) | day_night_1TRUE | day_night_2TRUE | day_night_3TRUE | day_night_4TRUE | day_night_5TRUE
      |   |   |   |   |                      1   | -3  | -0.43           | -0.64           | 1.5             | -2.3            | 2.6             ])
      |   |   |   |   '-- (1.28,7.54] (0.12 [ (0) | (I) 
      |   |   |   |                           1   | -2.3 ])
      |   |   |   '-- (1.28,7.54] (0.9 [ (0) | (I) 
      |   |   |                          1   | -2.1 ])
      |   |   '-- (1.28,7.54] (0.036 [ (0) | (I)  | day_night_1TRUE
      |   |                            1   | -2.7 | 1.3             ])
      |   '-- (1.28,7.54] (0.15 [ (0) | (I) 
      |                           1   | -1.6 ])
      '-- (1.28,7.54]
          +-- (0,1.28] (0.75 [ (0) | (I)
          |                    1   | 1.6 ])
          '-- (1.28,7.54]
              +-- (0,1.28] (1 [ (0) | (I)
              |                 1   | 1.8 ])
              '-- (1.28,7.54] (collapsing: 0.0032)
                  +-- (0,1.28] (0.029 [ (0) | (I)  | day_night_1TRUE
                  |                     1   | 0.99 | 1.3             ])
                  '-- (1.28,7.54] (0.72 [ (0) | (I)
                                          1   | 2.8 ])

---

    Code
      draw(t_m_cov_model, model = "coef", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ (0.48,1.63]/(0,0.48] | -2.7 | 1.1 
      |                       (1.63,7.54]/(0,0.48] | -4.5 | 0.84 ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ (0.48,1.63]/(0,0.48] | 0.12  | -9.1 | 9.6
      |   |                     (1.63,7.54]/(0,0.48] | -0.37 | -9.4 | 8.5 ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ (0.48,1.63]/(0,0.48] | 2     | 0.58
      |                                          (1.63,7.54]/(0,0.48] | -0.31 | 1    ])
      '-- (1.63,7.54] (2.4e-06 [ (0.48,1.63]/(0,0.48] | 10 | -7.5
                                 (1.63,7.54]/(0,0.48] | 13 | -8.4 ])

---

    Code
      draw(t_m_cov, model = "coef", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ (0.48,1.63]/(0,0.48] | -2.7 | 1.1 
      |                       (1.63,7.54]/(0,0.48] | -4.5 | 0.84 ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ (0.48,1.63]/(0,0.48] | 0.12  | -9.1 | 9.6
      |   |                     (1.63,7.54]/(0,0.48] | -0.37 | -9.4 | 8.5 ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ (0.48,1.63]/(0,0.48] | 2     | 0.58
      |                                          (1.63,7.54]/(0,0.48] | -0.31 | 1    ])
      '-- (1.63,7.54] (2.4e-06 [ (0.48,1.63]/(0,0.48] | 10 | -7.5
                                 (1.63,7.54]/(0,0.48] | 13 | -8.4 ])

---

    Code
      draw(t_m_cov_model, model = "full", time_sep = " | ", with_state = TRUE,
        digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ ((0,0.48])  | (I)  | day_night_1TRUE
      |                       (0.48,1.63] | -2.7 | 1.1            
      |                       (1.63,7.54] | -4.5 | 0.84            ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ ((0,0.48])  | (I)   | day_night_1TRUE | day_night_2TRUE
      |   |                     (0.48,1.63] | 0.12  | -9.1            | 9.6            
      |   |                     (1.63,7.54] | -0.37 | -9.4            | 8.5             ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ ((0,0.48])  | (I)   | day_night_1TRUE
      |                                          (0.48,1.63] | 2     | 0.58           
      |                                          (1.63,7.54] | -0.31 | 1               ])
      '-- (1.63,7.54] (2.4e-06 [ ((0,0.48])  | (I) | day_night_1TRUE
                                 (0.48,1.63] | 10  | -7.5           
                                 (1.63,7.54] | 13  | -8.4            ])

---

    Code
      draw(t_m_cov, model = "full", time_sep = " | ", with_state = TRUE, digits = 2)
    Output
      *
      +-- (0,0.48] (8.7e-07 [ ((0,0.48])  | (I)  | day_night_1TRUE
      |                       (0.48,1.63] | -2.7 | 1.1            
      |                       (1.63,7.54] | -4.5 | 0.84            ])
      '-- (0.48,1.63]
      |   +-- (0,0.48] (0.016 [ ((0,0.48])  | (I)   | day_night_1TRUE | day_night_2TRUE
      |   |                     (0.48,1.63] | 0.12  | -9.1            | 9.6            
      |   |                     (1.63,7.54] | -0.37 | -9.4            | 8.5             ])
      |   '-- (0.48,1.63], (1.63,7.54] (0.0049 [ ((0,0.48])  | (I)   | day_night_1TRUE
      |                                          (0.48,1.63] | 2     | 0.58           
      |                                          (1.63,7.54] | -0.31 | 1               ])
      '-- (1.63,7.54] (2.4e-06 [ ((0,0.48])  | (I) | day_night_1TRUE
                                 (0.48,1.63] | 10  | -7.5           
                                 (1.63,7.54] | 13  | -8.4            ])

