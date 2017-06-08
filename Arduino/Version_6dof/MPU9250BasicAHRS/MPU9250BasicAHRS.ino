/* MPU9250 Basic Example Code
  by: Kris Winer
  date: April 1, 2014
  license: Beerware - Use this code however you'd like. If you
  find it useful you can buy me a beer some time.
  Modified by Nicolas Rabany June 08, 2016

  Demonstrate basic MPU-9250 functionality including parameterizing the register
  addresses, initializing the sensor, getting properly scaled accelerometer,
  gyroscope, and magnetometer data out. Added display functions to allow display
  to on breadboard monitor. Addition of 9 DoF sensor fusion using open source
  Madgwick and Mahony filter algorithms. Sketch runs on the 3.3 V 8 MHz Pro Mini
  and the Teensy 3.1.

  SDA and SCL should have external pull-up resistors (to 3.3V).
  10k resistors are on the EMSENSR-9250 breakout board.

  Hardware setup:
  MPU9250 Breakout --------- Arduino
  VDD ---------------------- 3.3V
  VDDI --------------------- 3.3V
  SDA ----------------------- A4
  SCL ----------------------- A5
  GND ---------------------- GND
*/

#include "MPU9250.h"
#include <SdFat.h>
 
#define SerialDebug false  // Set to true to get Serial output for debugging

// Pin definitions
int intPin = 12;  // These can be changed, 2 and 3 are the Arduinos ext int pins
int myLed  = 13;  // Set up pin 13 led for toggling

MPU9250 myIMU;
unsigned long myTime=0;

// Replace SS with the chip slect pin for your SD.
const uint8_t sdChipSelect = 10;
// SD file system.
SdFat sd;
// Log file.
SdFile myFile;
uint32_t n = 0;
uint32_t t = 0;


void setup()
{
  Wire.begin();
  // TWBR = 12;  // 400 kbit/sec I2C speed
  Serial.begin(9600);

  /*---------IMU initialization---------------------------------------------*/ 
  // Set up the interrupt pin, its set as active high, push-pull
  pinMode(intPin, INPUT);
  digitalWrite(intPin, LOW);
  pinMode(myLed, OUTPUT);
  digitalWrite(myLed, HIGH);

  // Read the WHO_AM_I register, this is a good test of communication
  byte c = myIMU.readByte(MPU9250_ADDRESS, WHO_AM_I_MPU9250);
  Serial.print("MPU9250 "); Serial.print("I AM "); Serial.print(c, HEX);
  Serial.print(" I should be "); Serial.println(0x71, HEX);

  if (c == 0x71) // WHO_AM_I should always be 0x68
  {
    Serial.println("MPU9250 is online...");

    // Start by performing self test and reporting values
    myIMU.MPU9250SelfTest(myIMU.SelfTest);
    Serial.print("x-axis self test: acceleration trim within : ");
    Serial.print(myIMU.SelfTest[0], 1); Serial.println("% of factory value");
    Serial.print("y-axis self test: acceleration trim within : ");
    Serial.print(myIMU.SelfTest[1], 1); Serial.println("% of factory value");
    Serial.print("z-axis self test: acceleration trim within : ");
    Serial.print(myIMU.SelfTest[2], 1); Serial.println("% of factory value");
    Serial.print("x-axis self test: gyration trim within : ");
    Serial.print(myIMU.SelfTest[3], 1); Serial.println("% of factory value");
    Serial.print("y-axis self test: gyration trim within : ");
    Serial.print(myIMU.SelfTest[4], 1); Serial.println("% of factory value");
    Serial.print("z-axis self test: gyration trim within : ");
    Serial.print(myIMU.SelfTest[5], 1); Serial.println("% of factory value");


    myIMU.initMPU9250();
    // Initialize device for active mode read of acclerometer, gyroscope, and
    // temperature
    Serial.println("MPU9250 initialized for active data mode....");

    
  } // if (c == 0x71)
  else
  {
    Serial.print("Could not connect to MPU9250: 0x");
    Serial.println(c, HEX);
    while (1) ; // Loop forever if communication doesn't happen
  }
  /*---------end----------------------------------------------------*/ 

  /*---------microSD card initialization----------------------------*/ 
  Serial.begin(9600);
  // Initialize the SD and create or open the data file for append.
  if (!sd.begin(sdChipSelect, SPI_FULL_SPEED)
    || !myFile.open("data.txt", O_CREAT | O_WRITE | O_APPEND)) {
    // Replace this with somthing for your app.
    Serial.println(F("SD problem"));
    while(1);
  }
  myFile.println(" ");
  myFile.println("Initialization Success...");

  /*---------end----------------------------------------------------*/ 
}

void loop()
{
  delay(25);  // t(loop) ~= 25ms ==> delay of 25ms allows to recover data at freq of 20Hz
  // If intPin goes high, all data registers have new data
  // On interrupt, check if data ready interrupt
  if (myIMU.readByte(MPU9250_ADDRESS, INT_STATUS) & 0x01)
  {
    myIMU.readAccelData(myIMU.accelCount);  // Read the x/y/z adc values
    myIMU.getAres();

    // Get time of measurment
    myTime = millis();
    
    // Now we'll calculate the accleration value into actual g's
    // This depends on scale being set
    myIMU.ax = (float)myIMU.accelCount[0] * myIMU.aRes * (-1.0); // - accelBias[0];
    myIMU.ay = (float)myIMU.accelCount[1] * myIMU.aRes * (-1.0); // - accelBias[1];
    myIMU.az = (float)myIMU.accelCount[2] * myIMU.aRes * (-1.0); // - accelBias[2];

    // Convert in mg and add offset
    myIMU.ax = 1000 * myIMU.ax + 35.0;
    myIMU.ay = 1000 * myIMU.ay + 20.0;
    myIMU.az = 1000 * myIMU.az + 30.0;
    
    
    myIMU.readGyroData(myIMU.gyroCount);  // Read the x/y/z adc values
    myIMU.getGres();

    // Calculate the gyro value into actual degrees per second
    // This depends on scale being set
    myIMU.gx = (float)myIMU.gyroCount[0] * myIMU.gRes;
    myIMU.gy = (float)myIMU.gyroCount[1] * myIMU.gRes;
    myIMU.gz = (float)myIMU.gyroCount[2] * myIMU.gRes;

    // Add offset
    myIMU.gx = myIMU.gx + 0.300;
    myIMU.gy = myIMU.gy - 2.200;
    myIMU.gz = myIMU.gz - 0.600;
    
    
  } // if (readByte(MPU9250_ADDRESS, INT_STATUS) & 0x01)

 
  if (SerialDebug) //Writing in port monitor
  {
    // Print acceleration values in milligs!
    Serial.print("aX: "); Serial.print(myIMU.ax);
    Serial.print(" mg ");
    Serial.print("aY: "); Serial.print(myIMU.ay);
    Serial.print(" mg ");
    Serial.print("aZ: "); Serial.print(myIMU.az);
    Serial.println(" mg ");

    // Print gyro values in degree/sec
    Serial.print("gX: "); Serial.print(myIMU.gx, 3);
    Serial.print(" °/s ");
    Serial.print("gY: "); Serial.print(myIMU.gy, 3);
    Serial.print(" °/s ");
    Serial.print("gZ: "); Serial.print(myIMU.gz, 3);
    Serial.println(" °/s");
  } //(SerialDebug)
  
  else // Writting on microSD
  {
    myFile.print(myTime); // t
    myFile.print(" "); myFile.print(myIMU.ax); //aX (mg)
    myFile.print(" "); myFile.print(myIMU.ay); //aY (mg)
    myFile.print(" "); myFile.print(myIMU.az); //aZ (mg)
    myFile.print(" "); myFile.print(myIMU.gx, 3); //gX (°/sec)
    myFile.print(" "); myFile.print(myIMU.gy, 3); //gY (°/sec)
    myFile.print(" "); myFile.println(myIMU.gz, 3); //gZ (°/sec)
    
    // Use sync instead of close.
    myFile.sync();
  }
 
  
}

