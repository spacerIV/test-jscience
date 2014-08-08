import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import au.com.bytecode.opencsv.CSVReadProc;

import au.com.bytecode.opencsv.CSV;

import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileNotFoundException; 
import java.io.Writer;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.io.IOException;
import java.io.InputStreamReader;

import java.math.BigInteger;
import java.lang.Double;
import java.math.BigDecimal;
import java.util.Date;
import java.util.Random;
import java.util.Set;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;

import javax.measure.Measure;
import javax.measure.quantity.Mass;
import javax.measure.MeasureFormat;
import javax.measure.quantity.*;
import javax.measure.unit.*;

import org.jscience.economics.money.Currency;
import org.jscience.economics.money.Money;
import org.jscience.geography.coordinates.Altitude;
import org.jscience.geography.coordinates.CompoundCoordinates;
import org.jscience.geography.coordinates.LatLong;
import org.jscience.geography.coordinates.Time;
import org.jscience.geography.coordinates.UTM;
import org.jscience.geography.coordinates.XYZ;
import org.jscience.geography.coordinates.crs.CoordinatesConverter;
import org.jscience.mathematics.function.Polynomial;
import org.jscience.mathematics.function.Variable;
import org.jscience.mathematics.number.Complex;
import org.jscience.mathematics.number.Float64;
import org.jscience.mathematics.number.LargeInteger;
import org.jscience.mathematics.number.ModuloInteger;
import org.jscience.mathematics.number.Rational;
import org.jscience.mathematics.number.Real;
import org.jscience.mathematics.vector.ComplexMatrix;
import org.jscience.mathematics.vector.DenseMatrix;
import org.jscience.mathematics.vector.DenseVector;
import org.jscience.mathematics.vector.Float64Matrix;
import org.jscience.mathematics.vector.Matrix;
import org.jscience.mathematics.vector.Vector;
import org.jscience.physics.amount.Amount;
import org.jscience.physics.amount.AmountFormat;
import org.jscience.physics.model.RelativisticModel;

import javolution.lang.Configurable;
import javolution.lang.MathLib;
import javolution.text.TextBuilder;
import javolution.context.ConcurrentContext;
import javolution.context.LocalContext;
import javolution.context.StackContext;

import static javax.measure.unit.NonSI.*;
import static javax.measure.unit.SI.*;

import static org.jscience.economics.money.Currency.*;

public class Test implements javax.measure.quantity.Length
{

    public static void main(String[] args) throws IOException
    {
        Test test = new Test();

        //test.simple_convert();
        
        //test.getSIUnits();

        //test.getNonSIUnits();

        //test.test_meter();

        //test.test_volume();

        //test.getExact();

        //test.test_all();
        
        // test.more();

        //test.getConverter();

        //test.test_length();
        
        //test.test_length2();

        //test.test_length3();

        //test.test_format();

        //test.test_conversions();
        
        //test.test_acceleration();

        //test.test_hashmap();

        test.test_csv();

        //test.test_csv2();
    }

    public void test_conversions()
    {
        System.out.println( "" );
        System.out.println("Test Conversions");       

        Unit<Length> FURLONG = NonSI.FOOT.times(660);
        Unit<Length> MICRON = SI.MICRO(NonSI.INCH);

        //Derive the converter from the Unit
        System.out.println("Converter");       
        System.out.println(FURLONG.getConverterTo(SI.METER).convert(1));

        //Build a Scalar and have it calculate
        //System.out.println("Scalar");       
        //Quantity<Length> furlongScalar = new Scalar<Length>(1, FURLONG); // Ok.
        //final Unit<Length> mm = SI.MILLI(SI.METER);
        //System.out.println(furlongScalar.doubleValue(SI.METER));

        //Build a Measure and have it calculate
        System.out.println("Measure");       
        Measure lengthInFurlong = Measure.valueOf(1, FURLONG);
        System.out.println(lengthInFurlong.doubleValue(SI.METER));

        System.out.println(Measure.valueOf(1, FURLONG).doubleValue(SI.METER));
    }

    public void test_volume()
    {
        Amount<Volume> v1 = Amount.valueOf(20, 20, LITRE);
        Amount<Volume> v2 = Amount.valueOf(20, LITRE).to( MILLI( LITER ) );

        System.out.println( v1 );
        System.out.println( v2 );
    }

    public void test_meter()
    {
        System.out.println( "" );
        System.out.println("Test Meter");       

        //Amount<Length> mm = Amount.valueOf(1, METER).to( MILLI( METER ) );
        //Amount<Length> cm = Amount.valueOf(1, METER).to( CENTI( METER ) );
        //Amount<Length> pm = Amount.valueOf(1, METER).to( PICO( METER ) );
        //Amount<Length> yocto = Amount.valueOf( 1, METER).to( YOCTO( METER ) );

        Amount<Length> yotta = Amount.valueOf( 1, KILOMETER).to( YOTTA( METER ) );
        Amount<Length> zetta = Amount.valueOf( 1, KILOMETER).to( ZETTA( METER ) );
        Amount<Length> exa = Amount.valueOf( 1, KILOMETER).to( EXA( METER ) );
        Amount<Length> peta = Amount.valueOf( 1, KILOMETER).to( PETA( METER ) );
        Amount<Length> tera = Amount.valueOf( 1, KILOMETER).to( TERA( METER ) );
        Amount<Length> giga = Amount.valueOf( 1, KILOMETER).to( GIGA( METER ) );
        Amount<Length> mega = Amount.valueOf( 1, KILOMETER).to( MEGA( METER ) );
        Amount<Length> kilo = Amount.valueOf( 1, KILOMETER).to( KILO( METER ) );
        Amount<Length> hecto = Amount.valueOf( 1, KILOMETER).to( HECTO( METER ) );
        Amount<Length> deka = Amount.valueOf( 1, KILOMETER).to( DEKA( METER ) );
        Amount<Length> deci = Amount.valueOf( 1, KILOMETER).to( DECI( METER ) );
        Amount<Length> centi = Amount.valueOf( 1, KILOMETER).to( CENTI( METER ) );
        Amount<Length> milli = Amount.valueOf( 1, KILOMETER).to( MILLI( METER ) );
        Amount<Length> micro = Amount.valueOf( 1, KILOMETER).to( MICRO( METER ) );
        Amount<Length> nano = Amount.valueOf( 1, KILOMETER).to( NANO( METER ) );
        Amount<Length> pico = Amount.valueOf( 1, KILOMETER).to( PICO( METER ) );
        Amount<Length> femto = Amount.valueOf( 1, KILOMETER).to( FEMTO( METER ) );
        Amount<Length> atto = Amount.valueOf( 1, KILOMETER).to( ATTO( METER ) );
        Amount<Length> zepto = Amount.valueOf( 1, KILOMETER).to( ZEPTO( METER ) );
        Amount<Length> yocto = Amount.valueOf( 1, KILOMETER).to( YOCTO( METER ) );

        
        System.out.println("1 km = " + yotta );
        System.out.println("1 km = " + zetta );
        System.out.println("1 km = " + exa );
        System.out.println("1 km = " + peta );
        System.out.println("1 km = " + tera );
        System.out.println("1 km = " + giga );
        System.out.println("1 km = " + mega );
          System.out.println("---1 km = " + kilo );
          System.out.println("---1 km = " + hecto );
          System.out.println("---1 km = " + deka );
          System.out.println("---1 km = " + deci );
          System.out.println("---1 km = " + centi );
          System.out.println("---1 km = " + milli );
          System.out.println("---1 km = " + micro );
          System.out.println("---1 km = " + nano );
        System.out.println("1 km = " + pico );
        System.out.println("1 km = " + femto );
        System.out.println("1 km = " + atto );
        System.out.println("1 km = " + zepto );
        System.out.println("1 km = " + yocto ); 

    }

    public void getExact()
    {
        System.out.println("Exact Measurements");       
        Amount<Mass> m0 = Amount.valueOf(100, POUND);
        Amount<Mass> m1 = m0.times(33).divide(2);
        Amount<ElectricCurrent> m2 = Amount.valueOf("234 mA").to( MICRO(AMPERE) );
        
        System.out.println("m0 = " + m0);
        System.out.println("m1 = " + m1);
        System.out.println("m2 = " + m2);
    }

    public void getSIUnits() 
    {
        SI si = SI.getInstance();
        Set<Unit<?>> set_si_units = si.getUnits();

        System.out.println( "" );
        System.out.println( "SI units" );
        System.out.println( Arrays.toString( set_si_units.toArray() ) );
    }

    public void getNonSIUnits() 
    {
        NonSI nonSI = NonSI.getInstance();
        Set<Unit<?>> set_nonsi_units = nonSI.getUnits();

        System.out.println( "" );
        System.out.println( "NonSI units" );
        System.out.println( Arrays.toString( set_nonsi_units.toArray() ) );
    }

    public void simple_convert() 
    {
        final Unit<Angle> rev = NonSI.REVOLUTION.times(Math.PI * 2);
        final Unit<Length> cm = SI.CENTI(SI.METER);

        final Unit<? extends Quantity> cm3rev = cm.pow(3).divide(rev);
        final Unit<Length> mm = SI.MILLI(SI.METER);

        //double d1 = 4.44 cm3rev;
        Measure<Double, ? extends Quantity> d1 = Measure.valueOf(4.44, cm3rev);

        //double d2 = 2.22 mm;
        Measure<Double, Length> d2 = Measure.valueOf(2.22, mm);

        //double d = d1/d2;
        Unit<? extends Quantity> resultUnit = cm3rev.divide(mm);
        Measure<?, ?> result = Measure.valueOf(d1.getValue() / d2.getValue(), resultUnit);

        System.out.println("resultUnit "+resultUnit);
        System.out.println("result "+result);
    }

    public void test_all()
    {
        System.out.println("Testing...");   
        {
            System.out.println("");
            System.out.println("Exact Measurements");
            Amount<Mass> m0 = Amount.valueOf(100, POUND);
            Amount<Mass> m1 = m0.times(33).divide(2);
            Amount<ElectricCurrent> m2 = Amount.valueOf("234 mA").to(
                    MICRO(AMPERE));
            System.out.println("m0 = " + m0);
            System.out.println("m1 = " + m1);
            System.out.println("m2 = " + m2);

            System.out.println("");
            System.out.println("Inexact Measurements");
            Amount<Mass> m3 = Amount.valueOf(100.0, POUND);
            Amount<Mass> m4 = m0.divide(3);
            Amount<ElectricCurrent> m5 = Amount.valueOf("234 mA").to(AMPERE);
            Amount<Temperature> t0 = Amount.valueOf(-7.3, 0.5, CELSIUS);
            System.out.println("m3 = " + m3);
            System.out.println("m4 = " + m4);
            System.out.println("m5 = " + m5);
            System.out.println("t0 = " + t0);

            System.out.println("");
            System.out.println("Interval measurements");
            Amount<Volume> m6 = Amount.valueOf(20, 0.1, LITRE);
            Amount<Frequency> m7 = Amount.rangeOf(10, 11, KILO(HERTZ));
            System.out.println("m6 = " + m6);
            System.out.println("m7 = " + m7);

            System.out.println("");
            System.out.println("Amount.equals (identical) / Amount.approximates " +
                    "(takes into account errors such as numeric errors)");
            Amount<Frequency> m8 = Amount.valueOf(9000, HERTZ);
            Amount<Frequency> m10 = m8.divide(3).times(3); // Still exact.
            Amount<Frequency> m11 = m8.divide(7).times(7); // No more exact.
            System.out.println("m8 = " + m8);
            System.out.println("m10 = " + m10);
            System.out.println("m11 = " + m11);
            System.out.println("(m10 == m8) = " + m10.equals(m8));
            System.out.println("(m10 ≅ m8) = " + m10.approximates(m8));
            System.out.println("(m11 == m8) = " + m11.equals(m8));
            System.out.println("(m11 ≅ m8) = " + m11.approximates(m8));

            System.out.println("");
            System.out.println("AmountFormat - Plus/Minus Error (3 digits error)");
            AmountFormat.setInstance(AmountFormat
                    .getPlusMinusErrorInstance(3));
            System.out.println("m3 = " + m3);
            System.out.println("m4 = " + m4);
            System.out.println("m5 = " + m5);

            System.out.println("");
            System.out.println("AmountFormat - Bracket Error (2 digits error)");
            AmountFormat.setInstance(AmountFormat.getBracketErrorInstance(2));
            System.out.println("m3 = " + m3);
            System.out.println("m4 = " + m4);
            System.out.println("m5 = " + m5);

            System.out.println("");
            System.out.println("AmountFormat - Exact Digits Only");
            AmountFormat.setInstance(AmountFormat.getExactDigitsInstance());
            System.out.println("m3 = " + m3);
            System.out.println("m4 = " + m4);
            System.out.println("m5 = " + m5);

            System.out.println("");
            System.out.println("Numeric Errors");
            {
                Amount<Length> x = Amount.valueOf(1.0, METRE);
                Amount<Velocity> v = Amount.valueOf(0.01, METRES_PER_SECOND);
                Amount<Duration> t = Amount.valueOf(1.0, MICRO(SECOND));
                long ns = System.nanoTime();
                for (int i = 0; i < 10000000; i++) {
                    x = x.plus(v.times(t));
                }
                ns = System.nanoTime() - ns;
                AmountFormat.setInstance(AmountFormat
                        .getExactDigitsInstance());
                System.out.println(x
                        + " ("
                        + Amount.valueOf(ns, 0.5, NANO(SECOND)).to(
                                MILLI(SECOND)) + ")");
            }
            {
                double x = 1.0; // m
                double v = 0.01; // m/s
                double t = 1E-6; // s
                for (int i = 0; i < 10000000; i++) {
                    x += v * t; // Note: Most likely the compiler get v * t out of the loop.
                }
                System.out.println(x);
            }
            AmountFormat.setInstance(AmountFormat
                    .getPlusMinusErrorInstance(2));
        }
        {
            System.out.println("");
            System.out.println("Physical Models");
            // Selects a relativistic model for dimension checking (typically at start-up).
            RelativisticModel.select(); 

            // Length and Duration can be added.
            Amount<Length> x = Amount.valueOf(100, NonSI.INCH);
            x = x.plus(Amount.valueOf("2.3 µs")).to(METRE); 
            System.out.println(x); 
               
            // Energy is compatible with mass (E=mc2)
            Amount<Mass> m = Amount.valueOf("12 GeV").to(KILOGRAM); 
            System.out.println(m); 
        }

        {
            System.out.println("");
            System.out.println("Money/Currencies");
            ///////////////////////////////////////////////////////////////////////
            // Calculates the cost of a car trip in Europe for an American tourist.
            ///////////////////////////////////////////////////////////////////////

            // Use currency symbols instead of ISO-4217 codes.
            UnitFormat.getInstance().label(USD, "$"); // Use "$" symbol instead of currency code ("USD")
            UnitFormat.getInstance().label(EUR, "€"); // Use "€" symbol instead of currency code ("EUR")

            // Sets exchange rates.
            Currency.setReferenceCurrency(USD);
            EUR.setExchangeRate(1.17); // 1.0 € = 1.17 $

            // Calculates trip cost.
            Amount<?> carMileage = Amount.valueOf(20, MILE
                    .divide(GALLON_LIQUID_US)); // 20 mi/gal.
            Amount<?> gazPrice = Amount.valueOf(1.2, EUR.divide(LITRE)); // 1.2 €/L
            Amount<Length> tripDistance = Amount.valueOf(400, KILO(SI.METRE)); // 400 km
            Amount<Money> tripCost = tripDistance.divide(carMileage).times(
                    gazPrice).to(USD);
            // Displays cost.
            System.out.println("Trip cost = " + tripCost + " ("
                    + tripCost.to(EUR) + ")");
        }
        {
            System.out.println("");
            System.out.println("Matrices/Vectors");

            Amount<ElectricResistance> R1 = Amount.valueOf(100, 1, OHM); // 1% precision. 
            Amount<ElectricResistance> R2 = Amount.valueOf(300, 3, OHM); // 1% precision.
            Amount<ElectricPotential> U0 = Amount.valueOf(28, 0.01, VOLT); // ±0.01 V fluctuation.

            // Equations:  U0 = U1 + U2       |1  1  0 |   |U1|   |U0|
            //             U1 = R1 * I    =>  |-1 0  R1| * |U2| = |0 |
            //             U2 = R2 * I        |0 -1  R2|   |I |   |0 |
            //
            //                                    A      *  X   =  B
            //
            DenseMatrix<Amount<?>> A = DenseMatrix.valueOf(new Amount<?>[][] {
                { Amount.ONE,            Amount.ONE,            Amount.valueOf(0, OHM) },
                { Amount.ONE.opposite(), Amount.ZERO,           R1 },
                { Amount.ZERO,           Amount.ONE.opposite(), R2 } });
            DenseVector<Amount<?>> B = DenseVector.valueOf(new Amount<?>[] 
                { U0, Amount.valueOf(0, VOLT), Amount.valueOf(0, VOLT) });
            Vector<Amount<?>> X = A.solve(B);
            System.out.println(X);
            System.out.println(X.get(2).to(MILLI(AMPERE)));
        }
        {
            System.out.println("");
            System.out.println("Polynomials");

            // Defines two local variables (x, y).
            Variable<Complex> varX = new Variable.Local<Complex>("x");
            Variable<Complex> varY = new Variable.Local<Complex>("y");

            // f(x) = 1 + 2x + ix²
            Polynomial<Complex> x = Polynomial.valueOf(Complex.ONE, varX);
            Polynomial<Complex> fx = x.pow(2).times(Complex.I).plus(
                    x.times(Complex.valueOf(2, 0)).plus(Complex.ONE));
            System.out.println(fx);
            System.out.println(fx.pow(2));
            System.out.println(fx.differentiate(varX));
            System.out.println(fx.integrate(varY));
            System.out.println(fx.compose(fx));

            // Calculates expression.
            varX.set(Complex.valueOf(2, 3));
            System.out.println(fx.evaluate());
        }

        {
            System.out.println("");
            System.out.println("Coordinates Conversions");

            // Simple Lat/Long to UTM conversion.
            CoordinatesConverter<LatLong, UTM> latLongToUTM = LatLong.CRS
                    .getConverterTo(UTM.CRS);
            LatLong latLong = LatLong.valueOf(34.34, 23.56, DEGREE_ANGLE);
            UTM utm = latLongToUTM.convert(latLong);
            System.out.println(utm);

            // Lat/Long to XYZ conversion (assume height of zero).
            CoordinatesConverter<LatLong, XYZ> latLongToXYZ = LatLong.CRS
                    .getConverterTo(XYZ.CRS);
            XYZ xyz = latLongToXYZ.convert(latLong);
            System.out.println(xyz);

            // Compound coordinates - Lat/Long/Alt to XYZ conversion.
            Altitude alt = Altitude.valueOf(2000, FOOT);
            CompoundCoordinates<LatLong, Altitude> latLongAlt = 
                CompoundCoordinates.valueOf(latLong, alt);
            xyz = latLongAlt.getCoordinateReferenceSystem().getConverterTo(
                    XYZ.CRS).convert(latLongAlt);
            System.out.println(xyz);

            // Even more compounding...
            Time time = Time.valueOf(new Date());
            CompoundCoordinates<CompoundCoordinates<LatLong, Altitude>, Time> 
                latLongAltTime = CompoundCoordinates.valueOf(latLongAlt, time);
            System.out.println(latLongAltTime);
        }

        {
            System.out.println("");
            System.out.println("Numbers");

            Real two = Real.valueOf(2); // 2.0000..00 
            Real three = Real.valueOf(3);
            Real.setExactPrecision(100); // Assumes 100 exact digits for exact numbers.

            System.out.println("2/3       = " + two.divide(three));
            Real sqrt2 = two.sqrt();
            System.out.println("sqrt(2)   = " + sqrt2);
            System.out.println("Precision = " + sqrt2.getPrecision()
                    + " digits.");

            LargeInteger dividend = LargeInteger.valueOf("3133861182986538201");
            LargeInteger divisor = LargeInteger.valueOf("25147325102501733369");
            Rational rational = Rational.valueOf(dividend, divisor);
            System.out.println("rational  = " + rational);

            ModuloInteger m = ModuloInteger.valueOf("233424242346");
            LocalContext.enter(); // Avoids impacting others threads.
            try {
                ModuloInteger.setModulus(LargeInteger.valueOf("31225208137"));
                ModuloInteger inv = m.inverse();
                System.out.println("inverse modulo = " + inv);

                ModuloInteger one = inv.times(m);
                System.out.println("verification: one = " + one);

            } finally {
                LocalContext.exit();
            }

        }
    }

    public void more() 
    {
        System.out.println( "" );
        System.out.println( "More" );
        //Generate Units with scalar multiplication
        Unit<Length> FURLONG = NonSI.FOOT; // Ok
        System.out.println( FURLONG );
        FURLONG = NonSI.FOOT.times(2); // Ok
        System.out.println( FURLONG );
        // Unit<Length> FURLONG = SI.HOUR.times(60); // Compile error
        //Generate Units with prefixes
        Unit<Length> GIGAMETER = SI.GIGA(SI.METER); // Ok
        // Unit<Length> GIGAMETER = SI.GIGA(SI.NEWTON); // Compile error
        // Retrieval of the system unit (meter in this case)
        //System.out.println( NonSI.FOOT.getSystemUnit() );
        //System.out.println( NonSI.FOOT.getUnit() );
        // Retrieval of the unit dimension (L represents length)
        System.out.println(GIGAMETER.getDimension());
        // Dimension checking (allows/disallows conversions)
        System.out.println(NonSI.LIGHT_YEAR.isCompatible(NonSI.PARSEC)); //legal
        System.out.println(NonSI.LIGHT_YEAR.isCompatible(SI.WATT)); //illegal
        System.out.println(NonSI.LIGHT_YEAR.isCompatible(NonSI.PARSEC.divide(SI.SECOND))); //illegal
    }

    public void getConverter() 
    {
        System.out.println( "" );
        System.out.println( "getConverter" );

        Unit<Length> FURLONG = NonSI.FOOT.times(660);
        Unit<Length> MICROINCH = SI.MICRO(NonSI.INCH);
        System.out.println("FURLONG to MICROINCH: " + FURLONG.getConverterTo(MICROINCH).convert(1));
        System.out.println("FURLONG to METER:" +  SI.METER.getConverterTo(FURLONG).convert(1));

        System.out.println("MILE to METER: " + NonSI.MILE.getConverterTo(SI.METER).convert(1));

        System.out.println("MILE to KILOMETER: " + NonSI.MILE.getConverterTo(KILOMETER).convert(1));
    }

    public void test_format()
    {

        System.out.println("");
        System.out.println("Amount Format");
        Unit<Length> LEAGUE  = SI.KILOMETER.times(5.55600);
        Amount<Length> m = Amount.valueOf(1, METER).to( NonSI.LIGHT_YEAR ); 
        Amount<Mass> p = Amount.valueOf(100.0, POUND);
        //Measure<Double, Length> l = Measure.valueOf(100.0, SI.METER).to( NonSI.ANGSTROM );

        System.out.println("AmountFormat - Plus/Minus Error (3 digits error)");
        AmountFormat.setInstance(AmountFormat.getPlusMinusErrorInstance(3));
        //MeasureFormat.setInstance(MeasureFormat.getPlusMinusErrorInstance(3));
        //System.out.println( "meter to angstrom : " + l );
        //System.out.println( "meter to foot : "         + SI.METER.getConverterTo(  NonSI.FOOT ).convert( 1 ) );
        //System.out.println( "meter to league : "       + SI.METER.getConverterTo(  LEAGUE ).convert( 1 ) );
        System.out.println( "meter to lightyear : "       + m );
        System.out.println( "pound : "       + p );

        System.out.println("");
        System.out.println("AmountFormat - Bracket Error (2 digits error)");
        AmountFormat.setInstance(AmountFormat.getBracketErrorInstance(2));
        //MeasureFormat.setInstance(MeasureFormat.getPlusMinusErrorInstance(3));
        //System.out.println( "meter to angstrom : " + l );
        //System.out.println( "meter to foot : "         + SI.METER.getConverterTo(  NonSI.FOOT ).convert( 1 ) );
        //System.out.println( "meter to league : "       + SI.METER.getConverterTo(  LEAGUE ).convert( 1 ) );
        System.out.println( "meter to lightyear : "       + m );
        System.out.println( "pound : "       + p );

        System.out.println("");
        System.out.println("AmountFormat - Exact Digits Only");
        AmountFormat.setInstance(AmountFormat.getExactDigitsInstance());
        //MeasureFormat.setInstance(MeasureFormat.getPlusMinusErrorInstance(3));
        //System.out.println( "meter to angstrom : " + l );
        //System.out.println( "meter to foot : "         + SI.METER.getConverterTo(  NonSI.FOOT ).convert( 1 ) );
        //System.out.println( "meter to league : "       + SI.METER.getConverterTo(  LEAGUE ).convert( 1 ) );
        System.out.println( "meter to lightyear : "       + m );
        System.out.println( "pound : "       + p );

        Amount<Velocity> x = Amount.valueOf(7.5, NonSI.MILES_PER_HOUR);
        System.out.println(x);
        System.out.println( x.doubleValue(NonSI.MILES_PER_HOUR) + " miles per hour");
    }

    public void test_length() 
    {
        System.out.println( "" );
        System.out.println( "test length" );

        AmountFormat.setInstance(AmountFormat.getExactDigitsInstance());
        Unit<Length> FURLONG = NonSI.FOOT.times(660);
        Unit<Length> LEAGUE  = SI.METER.times(5556);

        System.out.println( "meter to angstrom : "     + SI.METER.getConverterTo(  NonSI.ANGSTROM ).convert( 1 ) );
        System.out.println( "meter to astronomical : " + SI.METER.getConverterTo(  NonSI.ASTRONOMICAL_UNIT ).convert( 1 ) );
        System.out.println( "meter to centimeter : "   + SI.METER.getConverterTo(  SI.CENTIMETER ).convert( 1 ) );
        //System.out.println( "meter to chain : "      + SI.METER.getConverterTo(  NonSI.CHAIN ).convert( 1 ) );
        //System.out.println( "meter to decimeter : "    + SI.METER.getConverterTo(  SI.DECIMETER ).convert( 1 ) );
        //System.out.println( "meter to fathom : "       + SI.METER.getConverterTo(  NonSI.FATHOM ).convert( 1 ) );
        System.out.println( "meter to foot : "         + SI.METER.getConverterTo(  NonSI.FOOT ).convert( 1 ) );
        System.out.println( "meter to furlong : "      + SI.METER.getConverterTo(  FURLONG ).convert( 1 ) );
        System.out.println( "meter to inch : "         + SI.METER.getConverterTo(  NonSI.INCH ).convert( 1 ) );
        System.out.println( "meter to kilometer : "    + SI.METER.getConverterTo(  SI.KILOMETER ).convert( 1 ) );
        System.out.println( "meter to league : "       + SI.METER.getConverterTo(  LEAGUE ).convert( 1 ) );
        System.out.println( "meter to light : "        + SI.METER.getConverterTo(  NonSI.LIGHT_YEAR ).convert( 1 ) );
        System.out.println( "meter to meter : "        + SI.METER.getConverterTo(  SI.METER ).convert( 1 ) );
        System.out.println( "meter to mile : "         + SI.METER.getConverterTo(  NonSI.MILE ).convert( 1 ) );
        System.out.println( "meter to millimeter : "   + SI.METER.getConverterTo(  SI.MILLIMETER ).convert( 1 ) );
        //System.out.println( "meter to micrometer : "   + SI.METER.getConverterTo(  SI.MICROMETER ).convert( 1 ) );
        //System.out.println( "meter to micron : "       + SI.METER.getConverterTo(  SI.MICRON ).convert( 1 ) );
        //System.out.println( "meter to nanometer : "    + SI.METER.getConverterTo(  SI.NANOMETER ).convert( 1 ) );
        System.out.println( "meter to nautical_mile: "     + SI.METER.getConverterTo(  NonSI.NAUTICAL_MILE ).convert( 1 ) );
        System.out.println( "meter to parsec : "       + SI.METER.getConverterTo(  NonSI.PARSEC ).convert( 1 ) );
        //System.out.println( "meter to rod : "          + SI.METER.getConverterTo(  NonSI.ROD ).convert( 1 ) );
        System.out.println( "meter to yard : "         + SI.METER.getConverterTo(  NonSI.YARD ).convert( 1 ) );

        //String str = new String("DEGREE_ANGLE");
        //UnitFormat u = u.getInstance();
        //System.out.println( " label: " + u.nameFor( SI.GRAM ) );
    }

    public void test_length2() 
    {
        System.out.println( "" );
        System.out.println( "test length2" );

        AmountFormat.setInstance(AmountFormat.getExactDigitsInstance());
        Unit<Length> FURLONG = NonSI.FOOT.times(660);
        Unit<Length> LEAGUE  = SI.METER.times(5556);
        
        System.out.println( "ANGSTROMS to METERS: "          + ( Measure.valueOf( 1, NonSI.ANGSTROM ).doubleValue( SI.METER) ) );
        System.out.println( "ASTRONOMICAL_UNITS to METERS: " + ( Measure.valueOf( 1, NonSI.ASTRONOMICAL_UNIT ).doubleValue( SI.METER) ) );
        System.out.println( "CENTIMETERS to METERS: "        + ( Measure.valueOf( 1, SI.CENTIMETER ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, NonSI.CHAIN ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, SI.DECIMETER ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, NonSI.FATHOM ).doubleValue( SI.METER) ) );
        System.out.println( "FOOTS to METERS: "              + ( Measure.valueOf( 1, NonSI.FOOT ).doubleValue( SI.METER) ) );
        System.out.println( "FURLONGS to METERS: "           + ( Measure.valueOf( 1, FURLONG ).doubleValue( SI.METER) ) );
        System.out.println( "INCHS to METERS: "              + ( Measure.valueOf( 1, NonSI.INCH ).doubleValue( SI.METER) ) );
        System.out.println( "KILOMETERS to METERS: "         + ( Measure.valueOf( 1, SI.KILOMETER ).doubleValue( SI.METER) ) );
        System.out.println( "LEAGUES to METERS: "            + ( Measure.valueOf( 1, LEAGUE ).doubleValue( SI.METER) ) );
        System.out.println( "LIGHT_YEARS to METERS: "        + ( Measure.valueOf( 1, NonSI.LIGHT_YEAR ).doubleValue( SI.METER) ) );
        System.out.println( "METERS to METERS: "             + ( Measure.valueOf( 1, SI.METER ).doubleValue( SI.METER) ) );
        System.out.println( "MILES to METERS: "              + ( Measure.valueOf( 1, NonSI.MILE ).doubleValue( SI.METER) ) );
        System.out.println( "MILLIMETERS to METERS: "        + ( Measure.valueOf( 1, SI.MILLIMETER ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, SI.MICROMETER ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, SI.MICRON ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, SI.NANOMETER ).doubleValue( SI.METER) ) );
        System.out.println( "NAUTICAL_MILES to METERS: "     + ( Measure.valueOf( 1, NonSI.NAUTICAL_MILE ).doubleValue( SI.METER) ) );
        System.out.println( "PARSECS to METERS: "            + ( Measure.valueOf( 1, NonSI.PARSEC ).doubleValue( SI.METER) ) );
        //System.out.println( ""                             + ( Measure.valueOf( 1, NonSI.ROD ).doubleValue( SI.METER) ) );
        System.out.println( "YARDS to METERS: "              + ( Measure.valueOf( 1, NonSI.YARD ).doubleValue( SI.METER) ) );
    }

    public void test_length3()
    {
        System.out.println( "" );
        System.out.println( "test length3" );

        System.out.println( "meter to mile1 : "         + SI.METER.getConverterTo(  NonSI.MILE ).convert( 1 ) );
        System.out.println( "meter to mile2 : "         + ( Measure.valueOf( 1, SI.METER ).doubleValue( NonSI.MILE ) ) );

        System.out.println( "mile to meter1 : "         + NonSI.MILE.getConverterTo(  SI.METER ).convert( 1 ) );
        System.out.println( "mile to meter2 : "         + ( Measure.valueOf( 1, NonSI.MILE ).doubleValue( SI.METER ) ) );

        System.out.println( "foot to mile1 : "         + NonSI.FOOT.getConverterTo(  NonSI.MILE ).convert( 1 ) );
        System.out.println( "foot to mile2 : "         + ( Measure.valueOf( 1, NonSI.FOOT ).doubleValue( NonSI.MILE ) ) );

        System.out.println( "mile to foot1 : "         + NonSI.MILE.getConverterTo(  NonSI.FOOT ).convert( 1 ) );
        System.out.println( "mile to foot2 : "         + ( Measure.valueOf( 1, NonSI.MILE ).doubleValue( NonSI.FOOT ) ) );

        // meter to mile1 : 6.213711922373339E-4
        // meter to mile2 : 6.213711922373339E-4
        // mile to meter1 : 1609.344
        // mile to meter2 : 1609.344
        //
        //convertpad
        //1 meter = 0.00062137 miles
        //1 mile  = 1609.34 meters 
       
        double dennis = 6.213711922373339E-4;
        System.out.println(dennis);
        System.out.println(String.format("%.7f", dennis));
        System.out.println(String.format("%.9f", new BigDecimal(dennis)));
        System.out.println(String.format("%.19f", new BigDecimal(dennis)));
        
    }

    public void test_acceleration() 
    {
        System.out.println("");
        System.out.println("test acceleration");

        Unit<Acceleration> CENTIMETERS_PER_SQUARE_SECOND = CENTI( SI.METERS_PER_SQUARE_SECOND );
        Unit<Acceleration> FOOT_PER_SQUARE_SECOND        = SI.METERS_PER_SQUARE_SECOND.times(0.30480);
        Unit<Acceleration> STANDARD_GRAVITY              = SI.METERS_PER_SQUARE_SECOND.times(9.80665);
        Unit<Acceleration> GAL                           = CENTIMETERS_PER_SQUARE_SECOND;
        Unit<Acceleration> MILLIGAL                      = MILLI( GAL );
        Unit<Acceleration> MICROGAL                      = MICRO( GAL );
        Unit<Acceleration> G_UNIT                        = STANDARD_GRAVITY;
        Unit<Acceleration> KILOMETERS_PER_SQUARE_HOUR    = SI.METERS_PER_SQUARE_SECOND.times(7.716049382716E-5);
        Unit<Acceleration> KILOMETERS_PER_HOUR_SECOND    = SI.METERS_PER_SQUARE_SECOND.times(0.2777777777778);
        Unit<Acceleration> MILES_PER_SQUARE_SECOND       = SI.METERS_PER_SQUARE_SECOND.times(1609.34348501);
        Unit<Acceleration> INCH_PER_SQUARE_SECOND        = SI.METERS_PER_SQUARE_SECOND.times(0.0254);
        Unit<Acceleration> YARD_PER_SQUARE_SECOND        = SI.METERS_PER_SQUARE_SECOND.times(0.9144);

        System.out.println("METERS_PER_SQUARE_SECOND:");
        System.out.println("1 METERS_PER_SQUARE_SECOND to CENTIMETERS_PER_SQUARE_SECOND: " + ( Measure.valueOf( 1, SI.METERS_PER_SQUARE_SECOND ).doubleValue( CENTI( SI.METERS_PER_SQUARE_SECOND)) ) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to FOOT_PER_SQUARE_SECOND: "        + SI.METERS_PER_SQUARE_SECOND.getConverterTo(FOOT_PER_SQUARE_SECOND).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to STANDARD_GRAVITY: "              + SI.METERS_PER_SQUARE_SECOND.getConverterTo(STANDARD_GRAVITY).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to GAL: "                           + SI.METERS_PER_SQUARE_SECOND.getConverterTo(GAL).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to G_UNIT: "                        + SI.METERS_PER_SQUARE_SECOND.getConverterTo(G_UNIT).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to KILOMETERS_PER_HOUR_SECOND: "    + SI.METERS_PER_SQUARE_SECOND.getConverterTo(KILOMETERS_PER_HOUR_SECOND).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to MILES_PER_SQUARE_SECOND: "       + SI.METERS_PER_SQUARE_SECOND.getConverterTo(MILES_PER_SQUARE_SECOND).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to INCH_PER_SQUARE_SECOND: "        + SI.METERS_PER_SQUARE_SECOND.getConverterTo(INCH_PER_SQUARE_SECOND).convert(1) );
        System.out.println("1 METERS_PER_SQUARE_SECOND to YARD_PER_SQUARE_SECOND: "        + SI.METERS_PER_SQUARE_SECOND.getConverterTo(YARD_PER_SQUARE_SECOND).convert(1) );

        System.out.println("FOOT_PER_SQUARE_SECOND:");
        System.out.println("1 FOOT_PER_SQUARE_SECOND to METERS_PER_SQUARE_SECOND: "      + FOOT_PER_SQUARE_SECOND.getConverterTo(METERS_PER_SQUARE_SECOND).convert(1) );
        System.out.println("1 FOOT_PER_SQUARE_SECOND to CENTIMETERS_PER_SQUARE_SECOND: " + FOOT_PER_SQUARE_SECOND.getConverterTo(CENTI(METERS_PER_SQUARE_SECOND)).convert(1) );
        System.out.println("1 INCH_PER_SQUARE_SECOND to METERS_PER_SQUARE_SECOND: "      + INCH_PER_SQUARE_SECOND.getConverterTo(METERS_PER_SQUARE_SECOND).convert(1) );
        System.out.println("1 YARD_PER_SQUARE_SECOND to METERS_PER_SQUARE_SECOND: "      + YARD_PER_SQUARE_SECOND.getConverterTo(METERS_PER_SQUARE_SECOND).convert(1) );
        System.out.println("1 YARD_PER_SQUARE_SECOND to INCH_PER_SQUARE_SECOND: "        + YARD_PER_SQUARE_SECOND.getConverterTo(INCH_PER_SQUARE_SECOND).convert(1) );
    }

    public void test_hashmap()
    {
        HashMap<String, String> symbols = new HashMap();

        symbols.put("G_UNIT", "g");
        symbols.put("GAL", "Gal");
        symbols.put("MILLIGAL", "mGal");
        symbols.put("MICROGAL", "\u00B5Gal");
        symbols.put("METERS_PER_SQUARE_SECOND", "m/s\u00B2");
        symbols.put("CENTIMETERS_PER_SQUARE_SECOND", "cm/s\u00B2");
        symbols.put("INCHES_PER_SQUARE_SECOND", "in/s\u00B2");
        symbols.put("FOOT_PER_SQUARE_SECOND", "ft/s\u00B2");
        symbols.put("YARDS_PER_SQUARE_SECOND", "yd/s\u00B2");
        symbols.put("MILES_PER_SQUARE_SECOND", "mi/s\u00B2");
        symbols.put("STANDARD_GRAVITY", "g\u2080");

        for ( String key : symbols.keySet() ) 
        {
            System.out.println( symbols.get(key) ); 
        }

        //System.out.println( symbols.get("METERS_PER_SQUARE_SECOND") );
    }

    public void test_csv() throws IOException
    {
        System.out.println("");
        System.out.println("test csv");

        CSVReader reader = null;
        String ADDRESS_FILE = "//Users/victor/Documents/github/test-jscience/res/acceleration.csv";

        try 
        {
            reader = new CSVReader(new FileReader(ADDRESS_FILE));
        }
        catch(  FileNotFoundException ex ) 
        { 
            System.out.println( "File not found exception, unfortunately" ); 
            throw(ex);
        } 

        TreeMap<String, String> unitSymbols             = new TreeMap<String, String>();
        TreeMap<String, String> unitNiceNames           = new TreeMap<String, String>();
        TreeMap<String, String> unitConversionUnit      = new TreeMap<String, String>();
        TreeMap<String, Double> unitTimes               = new TreeMap<String, Double>();
        TreeMap<String, Unit<Acceleration>> unitObjects = new TreeMap<String, Unit<Acceleration>>();

        String [] nextLine;
        while ((nextLine = reader.readNext()) != null) 
        {
            if ( nextLine[0].equals("name") && nextLine[1].equals("symbol") ) 
            {
                continue;
            }

            System.out.println("nicename: " + nextLine[2] );

            unitSymbols.put( nextLine[0], nextLine[1] );
            unitNiceNames.put( nextLine[0], nextLine[2] );
            unitConversionUnit.put( nextLine[0], nextLine[3] );
            unitTimes.put( nextLine[0], Double.parseDouble( nextLine[4] ) );
        }

        Unit<Acceleration> STANDARD_CONVERSION = SI.METERS_PER_SQUARE_SECOND;

        for ( String key : unitNiceNames.keySet() ) 
        {
            Unit<Acceleration> newObj = null; 
            if ( unitTimes.get( key ).equals( 1.0 ) ) 
            {
                newObj = STANDARD_CONVERSION;
            }
            else
            {
                newObj = STANDARD_CONVERSION.times( unitTimes.get(key) );
            }

            unitObjects.put( key, newObj );
        }

        System.out.println( "Done setup, now test: inch/sec\u00B2\u00B5\u2080 " );

        for ( String key : unitNiceNames.keySet() ) 
        {
            //System.out.println("1 INCH_PER_SQUARE_SECOND to METERS_PER_SQUARE_SECOND: "      + INCH_PER_SQUARE_SECOND.getConverterTo(METERS_PER_SQUARE_SECOND).convert(1) );
            Double d = unitObjects.get("CENTIMETERS_PER_SQUARE_SECOND").getConverterTo( unitObjects.get( key )).convert(1);  
            System.out.println( "1 " + unitNiceNames.get("CENTIMETERS_PER_SQUARE_SECOND") + " = " + d + " " + unitSymbols.get(key) );

        }
    }

    public void test_csv2() throws IOException
    {
        System.out.println("");
        System.out.println("test csv2");

        CSV csv = CSV.separator(',')
                     .quote('"')
                     .skipLines(1)
                     //.charset("UTF-8")
                     .create();

        String ADDRESS_FILE = "//Users/victor/Documents/github/test-jscience/res/acceleration.csv";

        TreeMap<String, String> unitSymbols             = new TreeMap<String, String>();
        TreeMap<String, String> unitNiceNames           = new TreeMap<String, String>();
        TreeMap<String, String> unitConversionUnit      = new TreeMap<String, String>();
        TreeMap<String, Double> unitTimes               = new TreeMap<String, Double>();
        TreeMap<String, Unit<Acceleration>> unitObjects = new TreeMap<String, Unit<Acceleration>>();

        csv.read(ADDRESS_FILE, new CSVReadProc() {
            public void procRow(int rowIndex, String... values) {
                if ( rowIndex == 1 )
                {
                    return;
                }
                System.out.println(rowIndex + ": " + Arrays.asList(values));
                System.out.println( values[0] );
                System.out.println( values[1] );
            }
        });
    }
}
