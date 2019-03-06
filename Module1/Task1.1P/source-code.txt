// Adapted from: https://learn.sparkfun.com/tutorials/force-sensitive-resistor-hookup-guide

const int FSR_PIN = A0; // Pin connected to FSR/resistor divider

const float VCC = 4.98; // Measured voltage of Ardunio 5V line
const float R_DIV = 3230.0; // Measured resistance of 3.3k resistor

void setup() 
{
  Serial.begin(9600);
  pinMode(FSR_PIN, INPUT);
  pinMode(13, OUTPUT);
}

void loop() 
{
  int fsrADC = analogRead(FSR_PIN);

  if (fsrADC != 0) // If the analog reading is non-zero
  {
    // Use ADC reading to calculate voltage:
    float fsrV = fsrADC * VCC / 1023.0;
    // Use voltage and static resistor value to 
    // calculate FSR resistance:
    float fsrR = R_DIV * (VCC / fsrV - 1.0);
    float force;
    float fsrG = 1.0 / fsrR; // Calculate conductance
    // Break parabolic curve down into two linear slopes:
    if (fsrR <= 600) 
      force = (fsrG - 0.00075) / 0.00000032639;
    else
      force =  fsrG / 0.000000642857;
    Serial.println("Force: " + String(force) + " g");
    Serial.println();

    if (force > 300)
      digitalWrite(13, HIGH);
    else
      digitalWrite(13, LOW);

    delay(500);
  }
  else
  {
    digitalWrite(13, LOW);
  }
}
