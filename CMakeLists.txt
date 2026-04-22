cmake_minimum_required (VERSION 3.10)
project ("C++_Linear_Algebra_LR1")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Собираем файлы напрямую из папки src
set(SOURCES
    src/main.cpp
    src/Matrix.cpp
    src/Solver.cpp
    src/Metrics.cpp
)

# Создаем экзешник в корне проекта
add_executable(C++_Linear_Algebra_LR1 ${SOURCES})

# Прописываем путь к хедерам
target_include_directories(C++_Linear_Algebra_LR1 PRIVATE src)