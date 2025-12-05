# Compilador y flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2 -Iinclude
LDFLAGS = 

# Directorios
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
RESULTS_DIR = results
SCRIPTS_DIR = scripts

# Archivos fuente y objeto
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
EXECUTABLE = $(BIN_DIR)/motor_bipedal

# Regla principal
all: crear_directorios $(EXECUTABLE)

# Crear directorios necesarios
crear_directorios:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR) $(RESULTS_DIR) $(SCRIPTS_DIR)

# Enlazar
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# Compilar
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ejecutar
run: all
	@echo "Ejecutando simulación..."
	@./$(EXECUTABLE)

# Visualizar (abrir gráficas)
view:
	@echo "Abriendo gráficas..."
	@if [ -f $(RESULTS_DIR)/posiciones.png ]; then \
		xdg-open $(RESULTS_DIR)/posiciones.png 2>/dev/null || open $(RESULTS_DIR)/posiciones.png 2>/dev/null || echo "  Abre manualmente: $(RESULTS_DIR)/posiciones.png"; \
	fi
	@sleep 1
	@if [ -f $(RESULTS_DIR)/velocidades.png ]; then \
		xdg-open $(RESULTS_DIR)/velocidades.png 2>/dev/null || open $(RESULTS_DIR)/velocidades.png 2>/dev/null || echo "  Abre manualmente: $(RESULTS_DIR)/velocidades.png"; \
	fi
	@sleep 1
	@if [ -f $(RESULTS_DIR)/animacion.gif ]; then \
		xdg-open $(RESULTS_DIR)/animacion.gif 2>/dev/null || open $(RESULTS_DIR)/animacion.gif 2>/dev/null || echo "  Abre manualmente: $(RESULTS_DIR)/animacion.gif"; \
	fi

# Limpiar
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(RESULTS_DIR)/* $(SCRIPTS_DIR)/*

# Limpiar solo resultados
clean_results:
	rm -f $(RESULTS_DIR)/* $(SCRIPTS_DIR)/*

.PHONY: all crear_directorios run view clean clean_results
