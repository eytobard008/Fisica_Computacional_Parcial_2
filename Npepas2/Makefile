CXX       = g++
CXXFLAGS  = -std=c++17 -Wall -O2
INCLUDES  = -Iinclude
SRC       = src/main.cpp src/Core.cpp src/Simulacion.cpp src/Visualizacion.cpp
TARGET    = bin/N_Particulas

$(TARGET): $(SRC)
	@mkdir -p bin
	@mkdir -p results
	@mkdir -p scripts
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

run: $(TARGET)
	./$(TARGET)

view:
	@echo "Abriendo visualizaciones..."
	@if [ -f results/Npepas.gif ]; then \
		echo "Abriendo animación principal..."; \
		(xdg-open results/Npepas.gif 2>/dev/null &) || (open results/Npepas.gif 2>/dev/null &) || echo "  Abre manualmente: results/Npepas.gif"; \
	else \
		echo "Error: results/Npepas.gif no existe."; \
	fi
	@sleep 1
	@if [ -f results/histograma.gif ]; then \
		echo "Abriendo histograma de velocidades..."; \
		(xdg-open results/histograma.gif 2>/dev/null &) || (open results/histograma.gif 2>/dev/null &) || echo "  Abre manualmente: results/histograma.gif"; \
	else \
		echo "Error: results/histograma.gif no existe."; \
	fi
	@sleep 1
	@if [ -f results/energia.png ]; then \
		echo "Abriendo gráfica de energía..."; \
		(xdg-open results/energia.png 2>/dev/null &) || (open results/energia.png 2>/dev/null &) || echo "  Abre manualmente: results/energia.png"; \
	else \
		echo "Error: results/energia.png no existe."; \
	fi
	@sleep 1
	@if [ -f results/presion.png ]; then \
		echo "Abriendo gráfica de presión..."; \
		(xdg-open results/presion.png 2>/dev/null &) || (open results/presion.png 2>/dev/null &) || echo "  Abre manualmente: results/presion.png"; \
	else \
		echo "Error: results/presion.png no existe."; \
	fi

clean:
	rm -f $(TARGET)
	rm -f results/*.dat results/*.gif results/*.png
	rm -f scripts/*.gnu

.PHONY: run clean view
