import java.util.PriorityQueue;
import java.util.Random;

/* Prints the reference output that each tests/jcompat/<name>_fixture.c must reproduce. */
public class JcompatFixtures {

    public static void main(String[] args) {
        switch (args[0]) {
            case "random" -> random();
            case "math" -> math();
            case "numbers" -> numbers();
            case "parse" -> parse();
            case "parseint" -> parseInt();
            case "utf8" -> utf8();
            case "pqueue" -> pqueue();
            case "search" -> search();
            default -> throw new IllegalArgumentException(args[0]);
        }
    }

    /* One line per evaluation: an operation, its input as raw bits, and the result. */
    static void numbers() {
        StringBuilder sb = new StringBuilder();
        Random r = new Random(20260925);
        java.util.List<Double> xs = new java.util.ArrayList<>(java.util.List.of(
                0.0, -0.0, 0.5, -0.5, 1.5, 2.5, -1.5, -2.5, 0.49999999999999994, -0.49999999999999994,
                4503599627370495.5, 4503599627370496.0, 2147483647.0, 2147483647.5, 2147483648.0,
                -2147483648.0, -2147483649.0, 9.223372036854775807E18, -9.223372036854776E18, 1e300, -1e300,
                Double.MIN_VALUE, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NaN,
                0.03125, 0.015625, 0.005, 0.015, 0.025, 0.125, 0.375, -0.001, -0.004, 0.995, 1.995, 0.00005));
        // Math.round edges and their double and float neighbours.
        double[] roundEdges = {0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.49999999999999994, 0x1p52, 0x1p52 + 1,
                -0x1p52, 0x1p53, 0x1p23, 0x1p23 + 1, 0x1p24, -0x1p23, 8388607.5, -8388607.5,
                4503599627370495.5, -4503599627370495.5, 0x1p31, -0x1p31, 0x1p63, -0x1p63, 0x1p64, -0x1p64,
                Float.MAX_VALUE, -Float.MAX_VALUE, Double.MAX_VALUE, -Double.MAX_VALUE};
        for (double e : roundEdges) {
            xs.addAll(java.util.List.of(e, Math.nextUp(e), Math.nextDown(e),
                    (double) Math.nextUp((float) e), (double) Math.nextDown((float) e)));
        }
        xs.addAll(java.util.List.of((double) 0.49999997f, (double) -0.49999997f, (double) Float.MIN_VALUE,
                (double) Float.NaN, (double) -0.0f, (double) Float.POSITIVE_INFINITY,
                (double) Float.NEGATIVE_INFINITY, (double) Integer.MAX_VALUE, (double) Integer.MIN_VALUE,
                (double) Long.MAX_VALUE, (double) Long.MIN_VALUE));
        for (int i = 0; i < 20000; ++i) {
            xs.add((double) r.nextFloat());
            xs.add((r.nextDouble() - 0.5) * Math.scalb(1.0, r.nextInt(70)));
        }
        for (double x : xs) {
            float fx = (float) x;
            sb.append("d2i ").append(d(x)).append(' ').append((int) x).append('\n');
            sb.append("d2l ").append(d(x)).append(' ').append((long) x).append('\n');
            sb.append("roundf ").append(f(fx)).append(' ').append(Math.round(fx)).append('\n');
            sb.append("roundd ").append(d(x)).append(' ').append(Math.round(x)).append('\n');
        }
        // DecimalFormat inputs in Beagle: j/100.0, and floats in [0, 1] widened to double.
        java.util.List<Double> formatted = new java.util.ArrayList<>();
        for (int j = 0; j <= 200; ++j) formatted.add(j/100.0);
        for (int k = 0; k <= 1 << 16; ++k) formatted.add((double) (k / 65536f));
        for (int i = 0; i < 20000; ++i) formatted.add((double) r.nextFloat());
        java.text.DecimalFormat hash2 = new java.text.DecimalFormat("#.##");
        java.text.DecimalFormat fixed2 = new java.text.DecimalFormat("0.00");
        java.text.DecimalFormat fixed4 = new java.text.DecimalFormat("0.0000");
        for (double x : formatted) {
            sb.append("hash2 ").append(d(x)).append(' ').append(hash2.format(x)).append('\n');
            sb.append("fixed2 ").append(d(x)).append(' ').append(fixed2.format(x)).append('\n');
            sb.append("fixed4 ").append(d(x)).append(' ').append(fixed4.format(x)).append('\n');
        }
        // Double.toString: specials, every power of two and ten and their
        // neighbours, the plain/scientific boundaries, and random values.
        java.util.List<Double> printed = new java.util.ArrayList<>(xs);
        printed.addAll(java.util.List.of(Double.MAX_VALUE, Double.MIN_NORMAL, 1e-3, 1e7, 1e23, 2e-3, 0.480887));
        for (int k = -1074; k <= 1023; ++k) printed.add(Math.scalb(1.0, k));
        for (int k = -324; k <= 308; ++k) printed.add(Double.parseDouble("1e" + k));
        for (int n = printed.size(), i = 0; i < n; ++i) {
            printed.add(Math.nextUp(printed.get(i)));
            printed.add(Math.nextDown(printed.get(i)));
        }
        for (int i = 0; i < 100000; ++i) {
            printed.add(Double.longBitsToDouble(r.nextLong()));
            printed.add(r.nextInt(100000000) / Math.pow(10, r.nextInt(12)));
        }
        for (double x : printed) {
            sb.append("tostring ").append(d(x)).append(' ').append(Double.toString(x)).append('\n');
        }
        // Float.toString: the same families in float, plus Par's bounds and the
        // values a user might type against them.
        java.util.List<Float> printedF = new java.util.ArrayList<>(java.util.List.of(
                0f, -0f, Float.MIN_VALUE, Float.MIN_NORMAL, Float.MAX_VALUE, Float.POSITIVE_INFINITY,
                Float.NEGATIVE_INFINITY, Float.NaN, 1e-3f, 1e7f, 0.5f, 0.6f, 1.5f, 2.5f, 0.1f, 0.005f,
                0.03125f, 0.480887f, 1.0e-45f, 3.4e38f, 9999999f, 1.0f, 100f, 1e6f, 2e-3f));
        for (int k = -149; k <= 127; ++k) printedF.add(Math.scalb(1.0f, k));
        for (int k = -45; k <= 38; ++k) printedF.add(Float.parseFloat("1e" + k));
        for (int n = printedF.size(), i = 0; i < n; ++i) {
            printedF.add(Math.nextUp(printedF.get(i)));
            printedF.add(Math.nextDown(printedF.get(i)));
        }
        for (int i = 0; i < 100000; ++i) {
            printedF.add(Float.intBitsToFloat(r.nextInt()));
            printedF.add((float) (r.nextInt(100000000) / Math.pow(10, r.nextInt(12))));
        }
        for (float x : printedF) {
            sb.append("ftostring ").append(f(x)).append(' ').append(Float.toString(x)).append('\n');
        }
        System.out.print(sb);
    }

    static String hex(byte[] b) {
        StringBuilder sb = new StringBuilder("x");
        for (byte x : b) sb.append(String.format("%02x", x & 0xff));
        return sb.toString();
    }

    /* One line per byte string without line terminators: its hex, then the UTF-8
     * re-encoding of Java's decoding by new String (BGZipIt) and by
     * InputStreamReader (BufferedReader.readLine). */
    static void utf8() {
        int[] pool = {0x41, 0x7f, 0x80, 0x8f, 0x90, 0x9f, 0xa0, 0xbf, 0xc0, 0xc1, 0xc2, 0xdf, 0xe0,
                0xe1, 0xec, 0xed, 0xee, 0xef, 0xf0, 0xf1, 0xf3, 0xf4, 0xf5, 0xf7, 0xf8, 0xff};
        Random r = new Random(20260926);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 200_000; ++i) {
            byte[] b = new byte[r.nextInt(9)];
            for (int j = 0; j < b.length; ++j) b[j] = (byte) pool[r.nextInt(pool.length)];
            byte[] viaString = new String(b, java.nio.charset.StandardCharsets.UTF_8)
                    .getBytes(java.nio.charset.StandardCharsets.UTF_8);
            String viaReader;
            try (java.io.Reader rd = new java.io.InputStreamReader(new java.io.ByteArrayInputStream(b),
                    java.nio.charset.StandardCharsets.UTF_8)) {
                StringBuilder s = new StringBuilder();
                for (int c; (c = rd.read()) != -1; ) s.append((char) c);
                viaReader = hex(s.toString().getBytes(java.nio.charset.StandardCharsets.UTF_8));
            } catch (java.io.IOException e) {
                throw new java.io.UncheckedIOException(e);
            }
            sb.append(hex(b)).append(' ').append(hex(viaString)).append(' ').append(viaReader).append('\n');
        }
        System.out.print(sb);
    }

    /* One line per string: its bytes in hex, then the bits of Double.parseDouble
     * and Float.parseFloat, or ERR where Java throws NumberFormatException. */
    static void parse() {
        Random r = new Random(20260927);
        java.util.List<String> texts = new java.util.ArrayList<>(java.util.List.of(
                "0", "-0", "+1", "1.", ".5", ".", "", " ", "1e", "1e+", "e5", "1e5", "1E-5", "1.5d", "1.5F", "1.5x",
                "NaN", "-NaN", "nan", "Infinity", "-Infinity", "+Infinity", "inf", "infinity", " 2.5 ", "\t3\n",
                "0x1p3", "0x1.8p1", "0X.8P0", "0x1", "0x1p", "0xp1", "0x1p3f", "1_000", "1,5", "--1", "+-1",
                "4.9e-324", "2.4703282292062328e-324", "2.4703282292062327e-324", "1.7976931348623157e308",
                "1.7976931348623159e308", "1e400", "1e-400", "3.4028235e38", "3.4028236e38", "1.4e-45", "7e-46",
                "0.1", "0.3", "123456789012345678901234567890", "9007199254740993", "1.00000017881393432617187499",
                "1.000000178813934326171875", "1.00000017881393432617187501"));
        String digits = "0123456789";
        for (int i = 0; i < 30000; ++i) {
            StringBuilder sb = new StringBuilder();
            if (r.nextInt(4) == 0) sb.append(r.nextBoolean() ? '-' : '+');
            int n = 1 + r.nextInt(25);
            for (int j = 0; j < n; ++j) sb.append(digits.charAt(r.nextInt(10)));
            if (r.nextBoolean()) {
                sb.insert(r.nextInt(sb.length() + 1), '.');
            }
            if (r.nextInt(3) == 0) sb.append(r.nextBoolean() ? 'e' : 'E').append(r.nextInt(700) - 350);
            if (r.nextInt(10) == 0) sb.append("fFdD".charAt(r.nextInt(4)));
            texts.add(sb.toString());
        }
        StringBuilder out = new StringBuilder();
        for (String t : texts) {
            StringBuilder hexText = new StringBuilder("x");
            for (byte b : t.getBytes(java.nio.charset.StandardCharsets.UTF_8)) hexText.append(String.format("%02x", b & 0xff));
            String dv, fv;
            try {
                dv = d(Double.parseDouble(t));
            } catch (NumberFormatException e) {
                dv = "ERR";
            }
            try {
                fv = f(Float.parseFloat(t));
            } catch (NumberFormatException e) {
                fv = "ERR";
            }
            out.append(hexText).append(' ').append(dv).append(' ').append(fv).append('\n');
        }
        System.out.print(out);
    }

    /* One line per byte string: its hex, then Integer.parseInt and
     * Long.parseLong of its UTF-8 decoding, each ERR where Java throws
     * NumberFormatException. */
    static void parseInt() {
        java.util.List<byte[]> inputs = new java.util.ArrayList<>();
        for (int c = 0; c <= 0xFFFF; ++c) {
            if (c >= 0xD800 && c <= 0xDFFF) continue;
            inputs.add(("7" + (char) c).getBytes(java.nio.charset.StandardCharsets.UTF_8));
        }
        for (String t : new String[] {"", "+", "-", "+0", "-0", "007", "2147483647", "2147483648",
                "-2147483648", "-2147483649", "99999999999", " 1", "1 ", "1_0", "0x10", "--1", "+-1",
                "٣٠٠٠", "-٣", "+٣", "１２３", "1٢३",
                "٩٩٩٩٩٩٩٩٩٩", "9223372036854775807", "9223372036854775808", "-9223372036854775808",
                "-9223372036854775809", "٩٢٢٣٣٧٢٠٣٦٨٥٤٧٧٥٨٠٧", "-٤٢"}) {
            inputs.add(t.getBytes(java.nio.charset.StandardCharsets.UTF_8));
        }
        for (int[] b : new int[][] {{0xf0, 0x9d, 0x9f, 0x8f}, {0x31, 0xf0, 0x9d, 0x9f, 0x8f}, {0xc0, 0xb1},
                {0xe0, 0x80, 0xb1}, {0xed, 0xa0, 0x80}, {0x80}, {0xd9}, {0x31, 0xd9},
                {0xf8, 0x88, 0x80, 0x80, 0x80}, {0xff}}) {
            byte[] bytes = new byte[b.length];
            for (int i = 0; i < b.length; ++i) bytes[i] = (byte) b[i];
            inputs.add(bytes);
        }
        StringBuilder out = new StringBuilder();
        for (byte[] bytes : inputs) {
            String t = new String(bytes, java.nio.charset.StandardCharsets.UTF_8);
            String v;
            try {
                v = Integer.toString(Integer.parseInt(t));
            } catch (NumberFormatException e) {
                v = "ERR";
            }
            String w;
            try {
                w = Long.toString(Long.parseLong(t));
            } catch (NumberFormatException e) {
                w = "ERR";
            }
            out.append(hex(bytes)).append(' ').append(v).append(' ').append(w).append('\n');
        }
        System.out.print(out);
    }

    static String d(double x) { return Long.toHexString(Double.doubleToRawLongBits(x)); }
    static String f(float x) { return Integer.toHexString(Float.floatToRawIntBits(x)); }

    /* One line per evaluation: an operation, its inputs as raw bits, and the result.
     * The strict-* lines check fdlibm itself. The site-* lines are Beagle's call
     * sites evaluated with Math, as Beagle does, on the platform running the test. */
    static void math() {
        StringBuilder sb = new StringBuilder();
        Random r = new Random(20260924);
        double[] specials = {0.0, -0.0, 1.0, -1.0, 2.0, 0.5, 10.0, 1e-310, Double.MIN_VALUE,
                Double.MAX_VALUE, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NaN};
        for (double x : specials) {
            sb.append("strict-log ").append(d(x)).append(' ').append(d(StrictMath.log(x))).append('\n');
            sb.append("strict-log10 ").append(d(x)).append(' ').append(d(StrictMath.log10(x))).append('\n');
            sb.append("strict-expm1 ").append(d(x)).append(' ').append(d(StrictMath.expm1(x))).append('\n');
            for (double y : specials) {
                sb.append("strict-pow ").append(d(x)).append(' ').append(d(y)).append(' ')
                        .append(d(StrictMath.pow(x, y))).append('\n');
            }
        }
        for (int i = 0; i < 20000; ++i) {
            double x = Math.scalb(r.nextDouble(), r.nextInt(200) - 100);
            double y = (r.nextDouble() - 0.5) * r.nextInt(64);
            double z = (r.nextDouble() - 0.5) * 80;
            sb.append("strict-log ").append(d(x)).append(' ').append(d(StrictMath.log(x))).append('\n');
            sb.append("strict-log10 ").append(d(x)).append(' ').append(d(StrictMath.log10(x))).append('\n');
            sb.append("strict-pow ").append(d(x)).append(' ').append(d(y)).append(' ')
                    .append(d(StrictMath.pow(x, y))).append('\n');
            sb.append("strict-expm1 ").append(d(z)).append(' ').append(d(StrictMath.expm1(z))).append('\n');
        }
        // main/Par.liStephensPMismatch
        for (int n = 1; n <= 200_000; n += 1 + n / 50) {
            double theta = 1/((Math.log(n) + 0.5));
            sb.append("site-pmismatch ").append(n).append(' ').append(f((float) (theta/(2*(theta + n))))).append('\n');
        }
        // bref/SeqCoder3.defaultMaxNSeq
        for (int n = 2; n <= 40_000_000; n += 1 + n / 100) {
            double exponent = 2*Math.log10(n) + 1;
            sb.append("site-maxnseq ").append(n).append(' ').append((long) Math.floor(Math.pow(2.0, exponent))).append('\n');
        }
        // phase/PhaseData.leaveUnphasedProp
        for (int it = 1; it <= 40; ++it) {
            for (int cnt = 0; cnt <= 5000; cnt += 1 + cnt / 20) {
                sb.append("site-unphased ").append(cnt).append(' ').append(it).append(' ')
                        .append(f((float) Math.pow(cnt, -1.0/it))).append('\n');
            }
        }
        // imp/ImpData.pRecomb and vcf/MarkerMap.pRecomb
        for (int i = 0; i < 20000; ++i) {
            double c = -(0.04*(float) (1 + r.nextInt(1_000_000))/(2 + r.nextInt(2_000_000)));
            double dist = Math.scalb(r.nextDouble(), -r.nextInt(30));
            sb.append("site-precomb ").append(d(c)).append(' ').append(d(dist)).append(' ')
                    .append(f((float) -Math.expm1(c*dist))).append('\n');
        }
        System.out.print(sb);
    }

    static final long[] SEEDS = {0L, 1L, -1L, -99999L, 42L, 0x5DEECE66DL,
            Long.MIN_VALUE, Long.MAX_VALUE, 1234567890123L, -1234567890123L};
    static final int[] BOUNDS = {1, 2, 3, 7, 10, 16, 181, 362, 1000, 1 << 20,
            (1 << 30) + 1, 1 << 30, Integer.MAX_VALUE, 1717986919};

    static void random() {
        StringBuilder sb = new StringBuilder();
        for (long seed : SEEDS) {
            Random r = new Random(seed);
            sb.append("seed ").append(seed).append('\n');
            for (int i = 0; i < 20; ++i) sb.append(r.nextInt()).append(' ');
            sb.append('\n');
            for (int bound : BOUNDS) {
                for (int i = 0; i < 20; ++i) sb.append(r.nextInt(bound)).append(' ');
                sb.append('\n');
            }
            for (int i = 0; i < 20; ++i) sb.append(r.nextLong()).append(' ');
            sb.append('\n');
            for (int i = 0; i < 64; ++i) sb.append(r.nextBoolean() ? '1' : '0');
            sb.append('\n');
            r.setSeed(seed + 7);
            for (int i = 0; i < 5; ++i) sb.append(r.nextInt(1000)).append(' ');
            sb.append('\n');
        }
        System.out.print(sb);
    }

    /* Ordered on key alone, like beagleutil.CompHapSegment on lastIbsStep, so
     * ties are broken by java.util.PriorityQueue's heap layout. */
    static final class Seg implements Comparable<Seg> {
        int key;
        final int id;
        Seg(int key, int id) {
            this.key = key;
            this.id = id;
        }
        @Override
        public int compareTo(Seg o) {
            return key != o.key ? (key < o.key ? -1 : 1) : 0;
        }
    }

    /* One line per Arrays.binarySearch(a, from, to, key): n, the n elements,
     * from, to, key and the result. Sorted arrays with runs of equal elements,
     * every subrange and every key from one below the minimum to one above the
     * maximum, so the probe order decides which equal element is found. */
    static void search() {
        StringBuilder out = new StringBuilder();
        Random r = new Random(7);
        java.util.List<int[]> arrays = new java.util.ArrayList<>();
        for (int n = 0; n <= 12; ++n) {
            for (int rep = 0; rep < 4; ++rep) {
                int[] a = new int[n];
                for (int i = 0; i < n; ++i) a[i] = r.nextInt(2 + rep * 3) - 2;
                java.util.Arrays.sort(a);
                arrays.add(a);
            }
        }
        arrays.add(new int[] {Integer.MIN_VALUE, Integer.MIN_VALUE, 0, Integer.MAX_VALUE, Integer.MAX_VALUE});
        for (int[] a : arrays) {
            java.util.TreeSet<Integer> keys = new java.util.TreeSet<>(java.util.List.of(0, -1, 1));
            for (int v : a) {
                keys.add(v);
                if (v > Integer.MIN_VALUE) keys.add(v - 1);
                if (v < Integer.MAX_VALUE) keys.add(v + 1);
            }
            for (int from = 0; from <= a.length; ++from) {
                for (int to = from; to <= a.length; ++to) {
                    for (int key : keys) {
                        out.append(a.length);
                        for (int v : a) out.append(' ').append(v);
                        out.append(' ').append(from).append(' ').append(to).append(' ').append(key).append(' ')
                                .append(java.util.Arrays.binarySearch(a, from, to, key)).append('\n');
                    }
                }
            }
        }
        System.out.print(out);
    }

    static void pqueue() {
        StringBuilder sb = new StringBuilder();
        for (long seed : SEEDS) {
            Random r = new Random(seed);
            PriorityQueue<Seg> q = new PriorityQueue<>(280);
            int nextId = 0;
            sb.append("seed ").append(seed).append('\n');
            for (int i = 0; i < 3000; ++i) {
                int op = r.nextInt(10);
                if (op < 5) {
                    q.offer(new Seg(r.nextInt(8), nextId++));
                } else if (op < 7) {
                    Seg s = q.poll();
                    sb.append(s == null ? -1 : s.id).append(' ');
                } else if (op == 7) {
                    Seg s = q.peek();
                    sb.append('p').append(s == null ? -1 : s.id).append(' ');
                } else if (op == 8) {
                    Seg s = q.poll();
                    if (s != null) {
                        s.key = r.nextInt(8);
                        q.offer(s);
                    }
                } else if (r.nextInt(50) == 0) {
                    q.clear();
                    sb.append("c ");
                }
            }
            sb.append('\n');
        }
        System.out.print(sb);
    }
}
