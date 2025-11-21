import java.util.Scanner;

public class GaussSeidelWithEquations {

    // ==========================
    // Clases y métodos auxiliares
    // ==========================

    private static class ParsedSide {
        double[] coeffs;   // coeficientes de las variables en este lado
        double constant;   // término independiente de este lado

        ParsedSide(double[] coeffs, double constant) {
            this.coeffs = coeffs;
            this.constant = constant;
        }
    }

    /**
     * Parsea una expresión lineal (un lado de la ecuación) y devuelve
     * los coeficientes de las variables y el término constante.
     *
     * Ejemplos de 'side':
     *  "4x-y+z"
     *  "2x+3y-5"
     *  "-x+2y+3"
     */
    private static ParsedSide parseSide(String side, String[] variables) {
        int n = variables.length;
        double[] coeffs = new double[n];
        double constant = 0.0;

        if (side == null || side.isEmpty()) {
            return new ParsedSide(coeffs, constant);
        }

        // Quitamos espacios y '*' por si escriben 2*x
        side = side.replace(" ", "").replace("*", "");

        if (side.isEmpty()) {
            return new ParsedSide(coeffs, constant);
        }

        // Aseguramos que empiece con signo
        if (side.charAt(0) != '+' && side.charAt(0) != '-') {
            side = "+" + side;
        }

        // Convertimos todos los '-' en '+-' para poder hacer split por '+'
        side = side.replace("-", "+-");

        // Separamos en términos
        String[] terms = side.split("\\+");

        for (String term : terms) {
            if (term == null || term.isEmpty()) {
                continue;
            }

            // ¿Este término tiene alguna variable?
            int varIndex = -1;
            String varNameFound = null;

            for (int i = 0; i < n; i++) {
                String var = variables[i];
                int idx = term.indexOf(var);
                if (idx != -1) {
                    varIndex = i;
                    varNameFound = var;
                    break;
                }
            }

            if (varIndex != -1) {
                // Término con variable, ej: "4x", "-y", "3.5z"
                String coefStr = term.replace(varNameFound, "");
                double coef;

                if (coefStr.isEmpty() || coefStr.equals("+")) {
                    coef = 1.0;
                } else if (coefStr.equals("-")) {
                    coef = -1.0;
                } else {
                    coef = Double.parseDouble(coefStr);
                }

                coeffs[varIndex] += coef;

            } else {
                // Término constante, ej: "5", "-3.2"
                double c = Double.parseDouble(term);
                constant += c;
            }
        }

        return new ParsedSide(coeffs, constant);
    }

    /**
     * Parsea una ecuación completa del tipo:
     *    "4x - y + z = 7"
     * usando los nombres de variables dados, y llena la fila 'row' de A y b.
     *
     * La ecuación general:
     *   LHS = RHS
     *
     * se reescribe como:
     *   (lhsCoeffs - rhsCoeffs) * x = rhsConst - lhsConst
     *
     * para obtener A[row][] y b[row].
     */
    private static void parseEquationToRow(String equation,
                                           String[] variables,
                                           double[][] A,
                                           double[] b,
                                           int row) {
        if (equation == null) {
            throw new IllegalArgumentException("La ecuación no puede ser nula.");
        }

        // Eliminamos espacios para simplificar
        equation = equation.replace(" ", "");

        // Separamos por '='
        String[] parts = equation.split("=");
        if (parts.length != 2) {
            throw new IllegalArgumentException(
                    "La ecuación debe tener un único signo '='. Ecuación: " + equation
            );
        }

        String left = parts[0];
        String right = parts[1];

        ParsedSide leftSide = parseSide(left, variables);
        ParsedSide rightSide = parseSide(right, variables);

        int n = variables.length;

        // A[row] = leftCoeffs - rightCoeffs
        for (int i = 0; i < n; i++) {
            A[row][i] = leftSide.coeffs[i] - rightSide.coeffs[i];
        }

        // b[row] = rightConst - leftConst
        b[row] = rightSide.constant - leftSide.constant;
    }

    // ==========================
    // Método de Gauss-Seidel
    // ==========================

    /**
     * Aplica el método de Gauss-Seidel para resolver A * x = b.
     *
     * @param A         Matriz de coeficientes (n x n)
     * @param b         Vector de términos independientes (n)
     * @param x0        Aproximación inicial (n)
     * @param tolerance Tolerancia para el criterio de parada
     * @param maxIter   Máximo número de iteraciones
     * @return Vector solución aproximada (n)
     */
    public static double[] gaussSeidel(double[][] A, double[] b, double[] x0,
                                       double tolerance, int maxIter) {
        int n = b.length;
        double[] x = new double[n];      // solución actual
        double[] xOld = new double[n];   // solución anterior

        // Copiamos x0 en x
        System.arraycopy(x0, 0, x, 0, n);

        for (int iter = 1; iter <= maxIter; iter++) {

            // Guardamos la solución anterior
            System.arraycopy(x, 0, xOld, 0, n);

            // Recorremos ecuación por ecuación
            for (int i = 0; i < n; i++) {
                double sum = 0.0;

                // sum = a[i][j] * x[j], j != i
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += A[i][j] * x[j];
                        // Para j < i usa x nuevo, para j > i usa x viejo (que aún está en x)
                    }
                }

                if (A[i][i] == 0.0) {
                    throw new ArithmeticException(
                            "Hay un cero en la diagonal en A[" + i + "][" + i + "]."
                    );
                }

                x[i] = (b[i] - sum) / A[i][i];
            }

            // Calculamos el error máximo entre x y xOld
            double maxError = 0.0;
            for (int i = 0; i < n; i++) {
                double error = Math.abs(x[i] - xOld[i]);
                if (error > maxError) {
                    maxError = error;
                }
            }

            System.out.printf("Iteración %d, error máximo = %.10f%n", iter, maxError);

            if (maxError < tolerance) {
                System.out.println("Convergencia alcanzada en " + iter + " iteraciones.");
                return x;
            }
        }

        System.out.println("No se alcanzó la convergencia en " + maxIter + " iteraciones.");
        return x;
    }

    /**
     * Comprueba si la matriz es diagonalmente dominante
     * (condición suficiente, pero no necesaria, para convergencia).
     */
    public static boolean isDiagonallyDominant(double[][] A) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            double diagonal = Math.abs(A[i][i]);
            double sumRow = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sumRow += Math.abs(A[i][j]);
                }
            }
            if (diagonal < sumRow) {
                return false;
            }
        }
        return true;
    }

    // ==========================
    // main: interacción con el usuario
    // ==========================

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.println("=== Método de Gauss-Seidel (entrada de ecuaciones) ===");
        System.out.print("Ingrese el número de ecuaciones/variables (n): ");

        // Leer n como línea completa para evitar problemas con saltos de línea
        int n = Integer.parseInt(scanner.nextLine().trim());

        System.out.println("\nIngrese los nombres de las variables en orden, separados por espacio.");
        System.out.println("Ejemplo para 3 variables: x y z");
        System.out.print("Variables: ");
        String[] variables = scanner.nextLine().trim().split("\\s+");

        if (variables.length != n) {
            throw new RuntimeException("numero de variables no coincide");
        }

        double[][] A = new double[n][n];
        double[] b = new double[n];

        System.out.println("\nIngrese cada ecuación en forma lineal, por ejemplo:");
        System.out.println("  4x - y + z = 7");
        System.out.println("  -2x + 6y + z = 9");
        System.out.println("  x + y + 5z = -6\n");

        for (int i = 0; i < n; i++) {
            System.out.print("Ecuación " + (i + 1) + ": ");
            String equation = scanner.nextLine();
            parseEquationToRow(equation, variables, A, b, i);
        }

        // Mostramos la matriz A y el vector b obtenidos
        System.out.println("\nMatriz A y vector b obtenidos a partir de las ecuaciones:");
        for (int i = 0; i < n; i++) {
            System.out.print("| ");
            for (int j = 0; j < n; j++) {
                System.out.printf("%10.4f ", A[i][j]);
            }
            System.out.printf("|   |x_%d|   =   %10.4f%n", i + 1, b[i]);
        }

        System.out.println();

        // Aproximación inicial: por defecto ceros
        double[] x0 = new double[n];
        System.out.println("Usaremos aproximación inicial x0 = (0, 0, ..., 0).");
        for (int i = 0; i < n; i++) {
            x0[i] = 0.0;
        }

        System.out.print("\nIngrese la tolerancia (ej. 1e-6): ");
        double tolerance = Double.parseDouble(scanner.nextLine().trim());

        System.out.print("Ingrese el máximo número de iteraciones: ");
        int maxIter = Integer.parseInt(scanner.nextLine().trim());

        // Comprobamos diagonalmente dominante
        if (!isDiagonallyDominant(A)) {
            System.out.println("\nADVERTENCIA: La matriz A no es diagonalmente dominante.");
            System.out.println("El método de Gauss-Seidel puede no converger.");
        }

        System.out.println("\nResolviendo el sistema con Gauss-Seidel...\n");
        double[] solution = gaussSeidel(A, b, x0, tolerance, maxIter);

        System.out.println("\nSolución aproximada:");
        for (int i = 0; i < solution.length; i++) {
            System.out.printf("%s = %.10f%n", variables[i], solution[i]);
        }

        scanner.close();
    }
}
