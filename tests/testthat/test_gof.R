# Create test data (this was simulated using dsims but for time and not having this as a dependency the data are pasted in here).

tot.data <- data.frame(Region.Label = "region",
                       Area = 1e+06,
                       Sample.Label = c(17, 17, 19, 19, 10, 10, 19, 19, 17, 17, 6, 6, 15, 15, 6, 6, 19, 19, 9, 9, 6, 6, 6, 6, 2, 2, 14, 14, 8, 8, 6, 6, 6, 6, 6, 6, 2, 2, 12, 12, 12, 12, 16, 16, 11, 11, 2, 2, 3, 3, 3, 3, 16, 16, 11, 11, 16, 16, 20, 20, 12, 12, 14, 14, 8, 8, 4, 4, 16, 16, 13, 13, 10, 10, 2, 2, 15, 15, 14, 14, 9, 9, 15, 15, 19, 19, 19, 19, 18, 18, 1, 1, 7, 7, 2, 2, 7, 7, 14, 14, 7, 7, 4, 4, 13, 13, 15, 15, 16, 16, 17, 17, 2, 2, 9, 9, 8, 8, 12, 12, 13, 13, 5, 5, 6, 6, 1, 1, 5, 5, 10, 10, 6, 6, 9, 9, 20, 20, 20, 20, 8, 8, 11, 11, 1, 1, 15, 15, 5, 5, 4, 4, 16, 16, 13, 13, 14, 14, 5, 5, 13, 13, 9, 9, 11, 11, 16, 16, 8, 8, 13, 13, 10, 10, 7, 7, 1, 1, 9, 9, 6, 6, 14, 14, 1, 1, 7, 7, 3, 3, 17, 17, 19, 19, 13, 13, 5, 5, 15, 15, 10, 10, 14, 14, 16, 16, 2, 2, 2, 2, 18, 18, 17, 17, 17, 17, 10, 10, 5, 5, 1, 1, 18, 18, 13, 13, 11, 11, 9, 9, 9, 9, 3, 3, 6, 6, 12, 12, 18, 18, 19, 19, 5, 5, 11, 11, 2, 2, 8, 8, 17, 17, 1, 1, 6, 6, 11, 11, 11, 11, 14, 14, 9, 9, 14, 14, 9, 9, 2, 2, 17, 17, 15, 15, 7, 7, 8, 8, 6, 6, 10, 10, 20, 20, 9, 9, 14, 14, 2, 2, 15, 15, 16, 16, 16, 16, 8, 8, 6, 6, 11, 11, 16, 16, 10, 10, 9, 9, 16, 16, 4, 4, 11, 11, 18, 18, 3, 3, 14, 14, 7, 7, 15, 15, 5, 5, 15, 15, 2, 2, 3, 3, 10, 10, 1, 1, 9, 9, 16, 16, 8, 8, 8, 8, 19, 19, 13, 13, 18, 18, 18, 18, 15, 15, 4, 4, 11, 11, 15, 15, 8, 8, 3, 3, 7, 7, 20, 20, 1, 1, 13, 13),
                       Effort = 500,
                       object = c(3, 3, 4, 4, 6, 6, 7, 7, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 33, 33, 34, 34, 35, 35, 36, 36, 40, 40, 41, 41, 43, 43, 46, 46, 47, 47, 48, 48, 50, 50, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 58, 58, 59, 59, 60, 60, 63, 63, 64, 64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 69, 70, 70, 71, 71, 72, 72, 74, 74, 75, 75, 76, 76, 77, 77, 78, 78, 81, 81, 82, 82, 83, 83, 84, 84, 85, 85, 86, 86, 87, 87, 89, 89, 90, 90, 92, 92, 93, 93, 94, 94, 96, 96, 97, 97, 98, 98, 99, 99, 100, 100, 102, 102, 103, 103, 104, 104, 106, 106, 107, 107, 110, 110, 111, 111, 112, 112, 113, 113, 114, 114, 115, 115, 117, 117, 118, 118, 119, 119, 120, 120, 121, 121, 122, 122, 124, 124, 125, 125, 128, 128, 129, 129, 131, 131, 132, 132, 133, 133, 134, 134, 135, 135, 136, 136, 137, 137, 139, 139, 140, 140, 141, 141, 142, 142, 144, 144, 145, 145, 146, 146, 148, 148, 149, 149, 150, 150, 154, 154, 155, 155, 156, 156, 157, 157, 158, 158, 159, 159, 160, 160, 161, 161, 162, 162, 163, 163, 164, 164, 165, 165, 166, 166, 168, 168, 169, 169, 170, 170, 172, 172, 173, 173, 175, 175, 176, 176, 177, 177, 178, 178, 179, 179, 180, 180, 182, 182, 184, 184, 185, 185, 187, 187, 188, 188, 190, 190, 192, 192, 193, 193, 196, 196, 197, 197, 198, 198, 201, 201, 202, 202, 203, 203, 204, 204, 205, 205, 206, 206, 207, 207, 208, 208, 209, 209, 210, 210, 212, 212, 213, 213, 216, 216, 217, 217, 218, 218, 219, 219, 220, 220, 222, 222, 223, 223, 224, 224, 225, 225, 226, 226, 227, 227, 229, 229, 230, 230, 231, 231, 232, 232, 234, 234, 236, 236, 240, 240, 241, 241, 242, 242, 245, 245, 246, 246, 247, 247, 249, 249),
                       observer = rep(c(1,2),186),
                       detected = c(1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                       distance = c(3.35, 3.35, 39.69, 39.69, 8.52, 8.52, 2.57, 2.57, 1.64, 1.64, 22.18, 22.18, 18.27, 18.27, 15.54, 15.54, 9.96, 9.96, 7.16, 7.16, 12.35, 12.35, 11.87, 11.87, 37.13, 37.13, 9.06, 9.06, 4.41, 4.41, 20.85, 20.85, 44.76, 44.76, 30, 30, 25.58, 25.58, 8.96, 8.96, 18.39, 18.39, 0.63, 0.63, 12.05, 12.05, 11.84, 11.84, 11.17, 11.17, 40.31, 40.31, 15.69, 15.69, 35.18, 35.18, 14.37, 14.37, 38.38, 38.38, 17.19, 17.19, 30.3, 30.3, 30.41, 30.41, 3.74, 3.74, 27.23, 27.23, 18.79, 18.79, 14.21, 14.21, 6, 6, 19.58, 19.58, 13.28, 13.28, 24.97, 24.97, 2.64, 2.64, 9.89, 9.89, 11.38, 11.38, 28.69, 28.69, 21.62, 21.62, 41.23, 41.23, 11.99, 11.99, 6.27, 6.27, 34.5, 34.5, 49.55, 49.55, 7.33, 7.33, 17.52, 17.52, 36.47, 36.47, 15.05, 15.05, 20.37, 20.37, 7.08, 7.08, 29.22, 29.22, 14.57, 14.57, 32.39, 32.39, 27, 27, 48.06, 48.06, 4.26, 4.26, 9.57, 9.57, 9.87, 9.87, 4.12, 4.12, 10.08, 10.08, 23.61, 23.61, 9.02, 9.02, 25.89, 25.89, 22.36, 22.36, 40.48, 40.48, 13.79, 13.79, 18.93, 18.93, 30, 30, 20.01, 20.01, 25.19, 25.19, 4.47, 4.47, 9.86, 9.86, 12.93, 12.93, 18.85, 18.85, 48.63, 48.63, 14.37, 14.37, 0.69, 0.69, 24.18, 24.18, 13.8, 13.8, 20.37, 20.37, 22.03, 22.03, 45.46, 45.46, 6.01, 6.01, 9.92, 9.92, 13.35, 13.35, 46.75, 46.75, 26.26, 26.26, 36.98, 36.98, 27.78, 27.78, 11.29, 11.29, 25.13, 25.13, 13.99, 13.99, 48.17, 48.17, 26.12, 26.12, 36.05, 36.05, 14.12, 14.12, 46.47, 46.47, 42.61, 42.61, 1.03, 1.03, 16.1, 16.1, 40.42, 40.42, 28.65, 28.65, 12.42, 12.42, 44.86, 44.86, 0.73, 0.73, 21.39, 21.39, 4.44, 4.44, 39.98, 39.98, 39.01, 39.01, 36.31, 36.31, 5.86, 5.86, 28.89, 28.89, 35.94, 35.94, 7.13, 7.13, 38.66, 38.66, 17.63, 17.63, 44.57, 44.57, 22.77, 22.77, 21.07, 21.07, 33.53, 33.53, 38.45, 38.45, 30.6, 30.6, 8.81, 8.81, 42.4, 42.4, 21.3, 21.3, 1.94, 1.94, 8.65, 8.65, 13.42, 13.42, 45.67, 45.67, 12.59, 12.59, 27.66, 27.66, 30.52, 30.52, 4.23, 4.23, 11.89, 11.89, 6.45, 6.45, 7.68, 7.68, 27.72, 27.72, 40.23, 40.23, 0.06, 0.06, 10.68, 10.68, 26.16, 26.16, 22.87, 22.87, 20, 20, 26.12, 26.12, 38.93, 38.93, 46.5, 46.5, 21.53, 21.53, 28.41, 28.41, 32.87, 32.87, 0.47, 0.47, 8.85, 8.85, 16.59, 16.59, 44.07, 44.07, 29, 29, 19.75, 19.75, 30.68, 30.68, 15.99, 15.99, 4.27, 4.27, 2.22, 2.22, 45.38, 45.38, 3.33, 3.33, 19.34, 19.34, 3.25, 3.25, 9.18, 9.18, 23.95, 23.95, 6.9, 6.9, 16.45, 16.45, 39.17, 39.17, 2.44, 2.44, 26.81, 26.81, 42.52, 42.52, 43.13, 43.13, 1.53, 1.53, 2.17, 2.17, 45.38, 45.38, 9.05, 9.05, 20.36, 20.36, 17.94, 17.94, 18.25, 18.25),
                       size = c(29, 29, 28, 28, 26, 26, 35, 35, 22, 22, 38, 38, 14, 14, 25, 25, 31, 31, 24, 24, 22, 22, 18, 18, 26, 26, 27, 27, 31, 31, 21, 21, 22, 22, 31, 31, 32, 32, 31, 31, 25, 25, 29, 29, 24, 24, 22, 22, 30, 30, 27, 27, 30, 30, 16, 16, 32, 32, 27, 27, 30, 30, 24, 24, 29, 29, 35, 35, 22, 22, 25, 25, 19, 19, 26, 26, 29, 29, 26, 26, 25, 25, 26, 26, 17, 17, 30, 30, 22, 22, 26, 26, 21, 21, 23, 23, 36, 36, 26, 26, 28, 28, 27, 27, 19, 19, 35, 35, 22, 22, 26, 26, 23, 23, 26, 26, 17, 17, 23, 23, 24, 24, 20, 20, 24, 24, 19, 19, 24, 24, 33, 33, 29, 29, 33, 33, 20, 20, 28, 28, 24, 24, 28, 28, 28, 28, 25, 25, 20, 20, 30, 30, 27, 27, 21, 21, 26, 26, 24, 24, 26, 26, 31, 31, 17, 17, 30, 30, 21, 21, 25, 25, 34, 34, 16, 16, 23, 23, 24, 24, 19, 19, 27, 27, 31, 31, 23, 23, 26, 26, 30, 30, 31, 31, 11, 11, 28, 28, 24, 24, 27, 27, 19, 19, 26, 26, 22, 22, 32, 32, 37, 37, 21, 21, 16, 16, 23, 23, 20, 20, 23, 23, 20, 20, 27, 27, 25, 25, 27, 27, 27, 27, 20, 20, 23, 23, 29, 29, 21, 21, 22, 22, 26, 26, 16, 16, 22, 22, 21, 21, 24, 24, 26, 26, 18, 18, 21, 21, 22, 22, 30, 30, 23, 23, 33, 33, 23, 23, 25, 25, 29, 29, 22, 22, 22, 22, 30, 30, 22, 22, 24, 24, 20, 20, 18, 18, 23, 23, 24, 24, 29, 29, 25, 25, 18, 18, 35, 35, 22, 22, 27, 27, 20, 20, 31, 31, 29, 29, 29, 29, 22, 22, 20, 20, 30, 30, 33, 33, 27, 27, 22, 22, 25, 25, 17, 17, 18, 18, 29, 29, 25, 25, 24, 24, 29, 29, 25, 25, 27, 27, 28, 28, 30, 30, 34, 34, 18, 18, 24, 24, 19, 19, 30, 30, 19, 19, 34, 34, 18, 18, 30, 30, 27, 27, 14, 14, 30, 30, 12, 12, 35, 35),
                       height = c(166, 166, 159, 159, 172, 172, 171, 171, 163, 163, 173, 173, 162, 162, 163, 163, 174, 174, 161, 161, 158, 158, 157, 157, 159, 159, 160, 160, 162, 162, 156, 156, 159, 159, 164, 164, 166, 166, 172, 172, 168, 168, 171, 171, 162, 162, 168, 168, 159, 159, 165, 165, 159, 159, 166, 166, 156, 156, 164, 164, 157, 157, 168, 168, 161, 161, 167, 167, 166, 166, 162, 162, 160, 160, 159, 159, 159, 159, 167, 167, 161, 161, 161, 161, 157, 157, 169, 169, 174, 174, 168, 168, 167, 167, 168, 168, 165, 165, 160, 160, 169, 169, 167, 167, 163, 163, 172, 172, 164, 164, 168, 168, 155, 155, 171, 171, 173, 173, 159, 159, 172, 172, 164, 164, 170, 170, 173, 173, 163, 163, 163, 163, 162, 162, 153, 153, 166, 166, 167, 167, 167, 167, 159, 159, 158, 158, 168, 168, 171, 171, 163, 163, 170, 170, 159, 159, 165, 165, 162, 162, 168, 168, 165, 165, 164, 164, 165, 165, 163, 163, 171, 171, 157, 157, 154, 154, 163, 163, 158, 158, 166, 166, 166, 166, 163, 163, 167, 167, 169, 169, 168, 168, 157, 157, 171, 171, 165, 165, 166, 166, 162, 162, 170, 170, 169, 169, 161, 161, 160, 160, 179, 179, 163, 163, 162, 162, 156, 156, 168, 168, 165, 165, 163, 163, 165, 165, 159, 159, 160, 160, 159, 159, 155, 155, 157, 157, 163, 163, 159, 159, 167, 167, 169, 169, 166, 166, 162, 162, 160, 160, 159, 159, 168, 168, 167, 167, 161, 161, 172, 172, 174, 174, 171, 171, 156, 156, 158, 158, 168, 168, 166, 166, 173, 173, 161, 161, 156, 156, 167, 167, 161, 161, 164, 164, 162, 162, 160, 160, 163, 163, 163, 163, 167, 167, 160, 160, 162, 162, 169, 169, 170, 170, 152, 152, 168, 168, 164, 164, 166, 166, 164, 164, 161, 161, 165, 165, 162, 162, 160, 160, 166, 166, 164, 164, 172, 172, 171, 171, 171, 171, 166, 166, 170, 170, 159, 159, 161, 161, 164, 164, 176, 176, 164, 164, 166, 166, 170, 170, 170, 170, 168, 168, 158, 158, 160, 160, 161, 161, 168, 168, 165, 165, 163, 163, 168, 168, 173, 173, 165, 165, 161, 161),
                       sex = c('male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'female', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'male', 'male', 'male'),
                       col = c('red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red'))

# Add distend distbegin to the data
tot.data <- Distance:::create_bins(tot.data, c(0,10,25,50))

# Code to simulate data above
# library(dsims)
# 
# covs <- list()
# covs$size  <- list(list(distribution = "poisson", lambda = 25))
# covs$height <- list(list(distribution = "normal", mean = 165, sd = 5))
# covs$sex   <- data.frame(level = c("male", "female"),
#                          prob = c(0.5, 0.5))
# covs$col   <- data.frame(level = c("blue", "red"),
#                          prob = c(0.4, 0.6))
# 
# # Define the population description (this time using the density to determine
# # the population size)
# popdesc <- make.population.description(covariates = covs,
#                                        N = 250)
# 
# cov.param <- list()
# cov.param$size <- c(log(1.02))
# cov.param$height <- c(log(1.005))
# cov.param$sex <- data.frame(level = c("male", "female"),
#                             param = c(log(1.5), 0))
# cov.param$col <- data.frame(level = c("blue", "red"),
#                             param = c(log(1.5), 0))
# 
# 
# # define the detecability
# detect <- make.detectability(key.function = "hn",
#                              scale.param = 5,
#                              cov.param = cov.param,
#                              truncation = 50)
# 
# plot(detect, popdesc)
# 
# simulation <- make.simulation(reps = 1,
#                               population.description = popdesc,
#                               detectability = detect)
# 
# # run an example survey to check the setup
# survey <- run.survey(simulation)
# 
# survey2 <- run.survey(survey, make.region())
# 
# ddata1 <- survey@dist.data
# ddata2 <- survey2@dist.data
# 
# head(ddata1)
# head(ddata2)
# 
# # Add observer and detected columns
# ddata1$observer <- 1
# ddata2$observer <- 2
# ddata1$detected <- 1
# ddata2$detected <- 1
# 
# # Combine the two datasets
# tot.data <- rbind(ddata1, ddata2)
# 
# # Check which objects only have 1 entry (all need 2)
# index <- which(table(tot.data$individual) == 1)
# index <- as.numeric(names(index))
# 
# # Add rows for animals that were detected by one observer but not the other
# tmp <- tot.data[tot.data$individual %in% index,]
# 
# # These animals were only detected by one observer
# tmp$detected <- 0
# new.observer <- tmp$observer
# new.observer <- ifelse(new.observer==1,2,1)
# tmp$observer <- new.observer
# 
# tot.data <- rbind(tot.data, tmp)
# 
# # Order the data - by observer
# index <- order(tot.data$observer)
# tot.data <- tot.data[index,]
# 
# # Order by individual
# index <- order(tot.data$individual)
# tot.data <- tot.data[index,]
# 
# # Remove object column
# tot.data <- tot.data[,-1]
# # Rename individual column to object
# names(tot.data)[1] <- "object"

test_that("gof when no df",{
  
  fit.mr <- ddf(data = tot.data,
                method = "io",
                dsmodel = ~mcds(key="hr", formula = ~size+height+sex),
                mrmodel = ~glm(formula = ~size+height+sex+col),
                meta.data = list(binned = TRUE, breaks = c(0,10,25,50)))
  
  # Check it knows there are insufficient df for pooled chisquare
  tmp <- ddf.gof(fit.mr)
  expect_true(is.na(tmp$chisquare$pooled.chi$df))
  expect_true(is.na(tmp$chisquare$pooled.chi$p))
  
})
