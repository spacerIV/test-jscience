
import java.math.BigInteger;
import java.util.Date;
import java.util.Random;

import javax.measure.Measure;
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


/**
import javax.measure.Measure;
import javax.measure.quantity.*;
import javax.measure.quantity.Angle;
import javax.measure.quantity.Length;
import javax.measure.quantity.Quantity;
import javax.measure.unit.NonSI;
import javax.measure.unit.SI;
import javax.measure.unit.Unit;

import javax.measure.quantity.*;
import javax.measure.unit.*;

import static javax.measure.unit.NonSI.*;
import static javax.measure.unit.SI.*;

import org.jscience.physics.amount.Amount;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
**/

public class Test
{
    public static void main(String[] args)
    {
        Test test = new Test();

        //test.simple_convert();
        
        //test.getNonSIUnits();

        //test.getExact();

        test.test_all();
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

    public void getNonSIUnits() 
    {
        //Set<Unit<?>> nonsi = SI.getUnits;
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
}
