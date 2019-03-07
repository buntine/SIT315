
const int FSR_PIN = 2;
const int LED_PIN = 13;

void pin_ISR() {
  int pinStatus = digitalRead(FSR_PIN);

  // I'm just tooling for any pressure here in order to enable the LED.
  if (pinStatus > 0) {
    digitalWrite(LED_PIN, HIGH);
    Serial.println("Turn on");
  } else {
    Serial.println("Turn off");
    digitalWrite(LED_PIN, LOW);
  }
}

void setup() 
{
  Serial.begin(9600);
  pinMode(FSR_PIN, INPUT);
  pinMode(13, OUTPUT);

  attachInterrupt(0, pin_ISR, CHANGE);
}

void loop() 
{
}
